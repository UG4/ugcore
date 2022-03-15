/*
 * Copyright (c) 2010-2022:  G-CSC, Goethe University Frankfurt
 * Author: Lukas Larisch
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__DISANCE_TO_BOUNDARY_BRUTEFORCE_2__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_SOLVER__DISANCE_TO_BOUNDARY_BRUTEFORCE_2__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/interface/matrix_operator_inverse.h"

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/function_spaces/dof_position_util.h"
#include "lib_disc/reference_element/reference_element.h"

#include "common/common.h"
#include "common/math/misc/math_util.h"

#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif
#include "../preconditioner/ilut_scalar.h"
#include "../interface/linear_operator_inverse.h"
#include "linear_solver.h"

#include "lib_algebra/cpu_algebra_types.h"

#include "../operator_util.h"

namespace ug{

template <typename TAlgebra, typename TDomain>
class DistanceToBoundaryBruteforce2
{
	public:
	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;
		typedef typename TAlgebra::matrix_type matrix_type;

		typedef MathVector<TDomain::dim> small_vector_type;

		/// Grid function type for the solution
		typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

		typedef typename GridFunc_t::template traits<Vertex>::const_iterator VertexConstIterator;

		/// Position attachment type
		typedef typename std::pair<small_vector_type, size_t> Position_t;

	public:
	///	constructor
		DistanceToBoundaryBruteforce2() : m_ssBoundaryIdx(-1), m_ssInnerIdx(-1){};


		const char* name() const {return "DistanceToBoundaryBruteforce2";}

	public:
		bool init(SmartPtr<GridFunc_t> spGF)
		{
			m_spGF = spGF;

			try{
				SmartPtr<DoFDistribution> dd = spGF->dof_distribution();

				m_vDistance = std::vector<number>(dd->num_indices(), 0.0);

				m_ssBoundaryIdx = spGF->domain()->subset_handler()->get_subset_index(m_ssBoundaryName);
				m_ssInnerIdx = spGF->domain()->subset_handler()->get_subset_index(m_ssInnerName);

				if(m_ssBoundaryIdx < 0)
					UG_THROW(name() << "::init_distance_to_boundary_bruteforce: Invalid boundary subset chosen! Use 'select_boundary(const char*)'.");

				if(m_ssInnerIdx < 0)
					UG_THROW(name() << "::init_distance_to_boundary_bruteforce: Invalid inner subset chosen! Use 'select_inner(const char*)'.");

				ConstSmartPtr<GridFunc_t> gf = spGF;
				ExtractPositions(gf->domain(), dd, m_vPositions);

			}UG_CATCH_THROW(name() << "::" << __FUNCTION__ << " failed")

			typedef typename GridFunc_t::template traits<ug::Vertex>::const_iterator ugVertIt_t;

			ugVertIt_t ugVertIt = spGF->template begin<ug::Vertex>();
			ugVertIt_t ugVertEnd = spGF->template end<ug::Vertex>();
			size_t k = 0;

			for(; ugVertIt != ugVertEnd; ++ugVertIt, ++k){
				ug::Vertex* v = *ugVertIt;

				int si = spGF->domain()->subset_handler()->get_subset_index(v);

				if(si == m_ssBoundaryIdx){
					m_vIdxBnd.push_back(k);
				}
				else if(si == m_ssInnerIdx){
					m_vIdxInner.push_back(k);
				}
				else{}
			}

			#ifdef UG_ENABLE_DEBUG_LOGS
			UG_LOG("Using " << name() << " (subset bnd " << m_ssBoundaryIdx << ", " << m_ssBoundaryName
						   << ", subset inner " << m_ssInnerIdx << ", " << m_ssInnerName 
						   << ", n = " << m_vIdxInner.size() << ", b = " << m_vIdxBnd.size() << ")\n");
			#endif

			return true;
		}

		bool solve()
		{
			try{

			
/*			//already done in init
			for(size_t j = 0; j < m_idxBnd.size(); ++j)
				m_vDistance[j] = 0.0f;
*/

			small_vector_type v;

			for(size_t i = 0; i < m_vIdxInner.size(); ++i){
				for(size_t j = 0; j < m_vIdxBnd.size(); ++j){
					VecSubtract(v, m_vPositions[m_vIdxInner[i]].first, m_vPositions[m_vIdxBnd[j]].first);
					number len = VecLength(v);
					m_vDistance[i] = (len < m_vDistance[i])? len : m_vDistance[i];
				}
			}

			std::vector<DoFIndex> ind(1);
			for(VertexConstIterator iter = m_spGF->template begin<Vertex>(); iter != m_spGF->template end<Vertex>(); ++iter)
			{
			//	get vertex
				Vertex* vrt = *iter;

			//	get vector holding all indices on the vertex
				m_spGF->inner_dof_indices(vrt, 0, ind);
				DoFRef(*m_spGF, ind[0]) = m_vDistance[ind[0][0]];
			}

			}UG_CATCH_THROW(name() << "::" << __FUNCTION__ << " failed")
			return true;
		}

		void select_boundary(const char* ssBoundaryName){
			m_ssBoundaryName = ssBoundaryName;
		}

		void select_inner(const char* ssInnerName){
			m_ssInnerName = ssInnerName;
		}

		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << name() << "\n";
			return ss.str();
		}

	///	Destructor
		virtual ~DistanceToBoundaryBruteforce2(){};

	protected:
		SmartPtr<GridFunc_t> m_spGF;
		int m_ssBoundaryIdx;
		int m_ssInnerIdx;
		const char* m_ssBoundaryName;
		const char* m_ssInnerName;

		std::vector<size_t> m_vIdxInner;
		std::vector<size_t> m_vIdxBnd;

		std::vector<number> m_vDistance;

		std::vector<Position_t> m_vPositions;

};

} // end namespace ug

#endif

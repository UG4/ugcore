/*
 * Copyright (c) 2011-2022:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__RIVERORDER__
#define __H__UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS__RIVERORDER__

#include <vector>
#include <utility> // for pair

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"
#include "lib_disc/function_spaces/dof_position_util.h"

#include "common/error.h"

namespace ug{

template <typename TAlgebra, typename TDomain, typename O_t>
class RiverOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	/// Grid function type for the solution
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;

	/// Position attachment type
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;

	RiverOrdering(){}

	/// clone constructor
	RiverOrdering( const LexOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new RiverOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	void compute(){
//topological sort here

		mat = NULL;
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	void init(M_t* A, const V_t& V){
//replace by subset check (not empty)
/*
		if(strcmp(m_dir, "") == 0){
			UG_THROW(name() << "::init: no direction chosen!");
		}
*/
		try{
			const GridFunc_t* pGridF;
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();

			size_t indN = dd->num_indices();

			if(indN != A->num_rows ()){
				UG_THROW(name() << "::init: #indices != #rows");
			}

			o.resize(indN);
			//ExtractPositions(pGridF->domain(), dd, m_vPositions);
		}
		catch(...){
			throw;
		}

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		mat = A;
	}

	void init(M_t*){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	virtual const char* name() const {return "RiverOrdering";}

	void select_sources(int subset){
		m_subset = subset;
	}

private:
	O_t o;
	M_t* mat;

	int m_subset;
	//std::vector<Position_t> m_vPositions;
};

} // end namespace ug

#endif

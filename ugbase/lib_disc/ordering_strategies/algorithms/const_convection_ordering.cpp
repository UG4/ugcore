/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_CONST_CONVECTION_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_CONST_CONVECTION_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/cuthill_mckee_ordering.hpp>

#include <set>
#include <algorithm> //reverse
#include <utility> //pair

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "bindings/lua/lua_user_data.h"

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"


namespace ug{

bool CompareScalar(const std::pair<number, size_t> &p1,
		const std::pair<number, size_t> &p2)
{
	return p1.first<p2.first;
}


template <typename TAlgebra, typename TDomain, typename O_t>
class ConstConvectionOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;
	typedef typename std::pair<number, size_t> Scalar_t;
	typedef MathVector<TDomain::dim> small_vec_t;
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;
	typedef UserData<MathVector<TDomain::dim>, TDomain::dim> Velocity_t;

	ConstConvectionOrdering(){}

	/// clone constructor
	ConstConvectionOrdering( const ConstConvectionOrdering<TAlgebra, TDomain, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new ConstConvectionOrdering<TAlgebra, TDomain, O_t>(*this));
	}

	void compute(){
		m_vScalars.resize(m_vPositions.size());
		small_vec_t pos;
		for(size_t i = 0; i < m_vPositions.size(); ++i){
			pos = m_vPositions[i].first;
			m_vScalars[i].first = pos*m_vel; //scalar product
			m_vScalars[i].second = m_vPositions[i].second;
		}

		std::sort(m_vScalars.begin(), m_vScalars.end(), CompareScalar);

		for(size_t i = 0; i < m_vScalars.size(); ++i){
			o[m_vScalars[i].second] = i;
		}

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void init(M_t* A, const V_t& V){
		const GridFunc_t* pGridF;
		unsigned n;

		try{
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();

			n = dd->num_indices();

			if(n != A->num_rows ()){
				UG_THROW(name() << "::init: #indices != #rows");
			}

			o.resize(n);
			ExtractPositions(pGridF->domain(), dd, m_vPositions);
		}
		catch(...){
			throw;
		}

		small_vec_t pos;
		pos = m_vPositions.at(0).first; //choose first occuring node

		//call lua function, store velocity in node with index 0 in m_vel
		(*m_spVelocity)(m_vel, pos, 0.0f, 0);

		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " (constant velocity " << m_vel << ")\n");
		#endif
	}

	void init(M_t* A){
		UG_THROW(name() << "::init: Cannot initialize smoother without a geometry. Specify the 2nd argument for init!");
	}

	void init(M_t* A, const V_t&, const O_t& inv_map){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t* A, const O_t& inv_map){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void check(){
		if(!is_permutation(o)){
			print(o);
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	void set_velocity(const char* strVelocity){
		m_spVelocity = make_sp(new LuaUserData<MathVector<TDomain::dim>, TDomain::dim>(strVelocity));
	}

	virtual const char* name() const {return "FollowConvectionOrdering";}

private:
	O_t o;

	SmartPtr<Velocity_t> m_spVelocity;
	small_vec_t m_vel;

	std::vector<Position_t> m_vPositions;
	std::vector<Scalar_t> m_vScalars;
};

} //namespace


#endif //guard

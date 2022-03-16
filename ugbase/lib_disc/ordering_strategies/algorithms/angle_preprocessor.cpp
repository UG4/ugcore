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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_ANGLE_PREPROCESSOR__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_ANGLE_PREPROCESSOR__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <vector>
#include <utility> //for pair
#include <climits> //INT_MAX

#include "lib_algebra/ordering_strategies/algorithms/IOrderingPreprocessor.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/grid_function.h"

#include "lib_disc/spatial_disc/user_data/user_data.h"
#include "bindings/lua/lua_user_data.h"

#include <boost/graph/connected_components.hpp>

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{

template <typename TAlgebra, typename TDomain, typename G_t>
class AnglePreprocessor : public IOrderingPreprocessor<TAlgebra, G_t>
{
public:
	typedef typename TAlgebra::vector_type V_t;
	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd_t;
	typedef typename boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator nIt_t;
	typedef std::pair<vd_t, int> indegs_t;
	typedef IOrderingPreprocessor<TAlgebra, G_t> baseclass;
	typedef typename std::pair<MathVector<TDomain::dim>, size_t> Position_t;
	typedef MathVector<TDomain::dim> small_vec_t;
	typedef GridFunction<TDomain, TAlgebra> GridFunc_t;
	typedef UserData<MathVector<TDomain::dim>, TDomain::dim> Velocity_t;

	AnglePreprocessor(){}

	/// clone constructor
	AnglePreprocessor( const AnglePreprocessor<TAlgebra, TDomain, G_t> &parent )
			: baseclass(), m_threshold(0.){}

	SmartPtr<IOrderingPreprocessor<TAlgebra, G_t> > clone()
	{
		return make_sp(new AnglePreprocessor<TAlgebra, TDomain, G_t>(*this));
	}

	void compute(){
		unsigned n = boost::num_vertices(*g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
		}

		if(m_threshold == 0){
			UG_THROW(name() << "::compute: No threshold set! Use 'set_threshold(number angle)'");
		}

		//eIt_t eIt, eEnd;
		vIt_t vIt, vEnd;
		vd_t s, t;
		small_vec_t pos_s, pos_t, vel_s, dir_st;
		nIt_t nIt, nEnd;
		size_t numRemoved = 0;
		number angle;


		std::vector<int> component(boost::num_vertices(*g));
		int numCompBegin = boost::connected_components(*g, &component[0]);
		

		for(boost::tie(vIt, vEnd) = boost::vertices(*g); vIt != vEnd; ++vIt){
			s = *vIt;
			pos_s = m_vPositions.at(s).first;

			//call lua function, store velocity of s in vel_s
			(*m_spVelocity)(vel_s, pos_s, 0.0f, 0);

			if(VecLengthSq(vel_s) > 0){
				std::vector<vd_t> to_delete;
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, *g); nIt != nEnd; ++nIt){
					t = *nIt;
					pos_t = m_vPositions.at(t).first;
					VecSubtract(dir_st, pos_s, pos_t);

					angle = VecAngle(dir_st, vel_s);
					if(angle > m_threshold){
						boost::remove_edge(s, t, *g);
						++numRemoved;
					}


/*
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, g); nIt != nEnd; ++nIt){
					t = *nIt;
					pos_t = m_vPositions.at(t).first;
					VecSubtract(dir_st, pos_s, pos_t);

					angle = VecAngle(dir_st, vel_s);
					if(angle > threshold){
						boost::remove_edge(s, t, g);
						++numRemoved;
					}
*/
				}
			}
			else{
				numRemoved += boost::degree(s, *g);
				boost::clear_vertex(*vIt, *g);
			}
		}

		int numCompEnd = boost::connected_components(*g, &component[0]);

		UG_LOG(name() << "::compute: " << numRemoved << " edges removed. Remaining: " << boost::num_edges(*g) << "\n");
		UG_LOG(name() << "::compute: " << "#components " << numCompBegin << " -> " << numCompEnd << "\n");
	}

	void init(const V_t &V, G_t &G){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		g = &G;

		try{
			const GridFunc_t* pGridF;
			if((pGridF = dynamic_cast<const GridFunc_t*>(&V)) == 0){
				UG_THROW(name() << "::init: No DoFDistribution specified.");
			}

			SmartPtr<DoFDistribution> dd = ((GridFunc_t*) pGridF)->dof_distribution();
			ExtractPositions(pGridF->domain(), dd, m_vPositions);
		}
		catch(...){
			throw;
		}
	}

	void set_velocity(const char* strVelocity){
		m_spVelocity = make_sp(new LuaUserData<MathVector<TDomain::dim>, TDomain::dim>(strVelocity));
	}

	void set_threshold(number threshold){
		m_threshold = threshold;
	}

	virtual const char* name() const {return "AnglePreprocessor";}

private:
	G_t *g;

	std::vector<size_t> indegs;

	std::vector<Position_t> m_vPositions;
	SmartPtr<Velocity_t> m_spVelocity;

	number m_threshold;
};

} //namespace

#endif //guard

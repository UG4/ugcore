/*
 * Copyright (c) 2021:  G-CSC, Goethe University Frankfurt
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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_TOPOLOGICAL_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_TOPOLOGICAL_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <vector>
#include <utility> //for pair
#include <climits> //INT_MAX

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

#include "lua_ordering.h"

namespace ug{

//for cycle-free matrices only
template <typename TAlgebra, typename O_t>
class TopologicalOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;
	typedef boost::graph_traits<G_t>::vertex_descriptor vd_t;
	typedef boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef boost::graph_traits<G_t>::adjacency_iterator nIt_t;
	typedef std::pair<vd_t, int> indegs_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	TopologicalOrdering(){}

	/// clone constructor
	TopologicalOrdering( const TopologicalOrdering<TAlgebra, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new TopologicalOrdering<TAlgebra, O_t>(*this));
	}

	indegs_t min_indegree_vertex(){
		vd_t minv;
		int mind = INT_MAX;
		vIt_t vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			int d = indegs[*vIt];
			if(d >= 0 && d < mind){
				mind = d;
				minv = *vIt;
			}
		}
		return std::make_pair(minv, mind);
	}

	//TODO: this is a very inefficient implementation, rewrite if necessary
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		o.resize(n);
		indegs.resize(n);

		vIt_t vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			indegs[*vIt] = boost::in_degree(*vIt, g);
		}

		nIt_t nIt, nEnd;
		for(unsigned i = 0; i < n; ++i){
			indegs_t indeg = min_indegree_vertex();
			if(indeg.second == 0){
				o[indeg.first] = i;
				--indegs[indeg.first]; //becomes -1
				for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(indeg.first, g); nIt != nEnd; ++nIt){
					--indegs[*nIt];
				}
			}
			else{
				UG_THROW(name() << "::compute: Graph is not cycle-free!");
			}
		}

		g = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void init(M_t* A, const V_t&){
		init(A);
	}

	void init(M_t* A){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		m_A = A;

		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					boost::add_edge(i, conn.index(), g);
				}
			}
		}
	}

	void init(M_t*, const V_t&, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	void init(M_t*, const O_t&){
		UG_THROW(name() << "::init: Algorithm does not support induced subgraph version!");
	}

	void check(){
		UG_COND_THROW(!is_permutation(o), name() << "::check: Not a permutation!");
	}

	O_t& ordering(){
		return o;
	}

	SmartPtr<LuaOrdering> get_lua_ordering(){
		SmartPtr<LuaOrdering> lua_ord = SmartPtr<LuaOrdering>(new LuaOrdering());
		lua_ord->ordering = o;
		return lua_ord;
	}

	virtual const char* name() const {return "TopologicalOrdering";}

private:
	G_t g;
	O_t o;
	M_t* m_A;

	std::vector<size_t> indegs;
};

} //namespace

#endif //guard

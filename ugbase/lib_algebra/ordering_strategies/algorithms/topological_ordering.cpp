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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_TOPOLOGICAL_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_TOPOLOGICAL_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <vector>
#include <utility> //for pair
#include <deque>

#include "IOrderingAlgorithm.h"
#include "util.cpp"

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{

template <typename O_t, typename G_t>
void topological_ordering_core(O_t& o, G_t& g, bool inverse){
	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd_t;
	typedef typename boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator nIt_t;

	size_t n = boost::num_vertices(g);

	if(n == 0){
		UG_THROW("topological_ordering_core: Graph is empty!");
	}

	o.resize(n);
	std::vector<int> indegs(n);
	std::vector<BOOL> isindeq(n, false);
	std::deque<vd_t> deq;

	vIt_t vIt, vEnd;
	nIt_t nIt, nEnd;
	for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
		for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(*vIt, g); nIt != nEnd; ++nIt){
			++indegs[*nIt];
		}
	}

	for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
		if(indegs[*vIt] == 0){
			deq.push_front(*vIt);
			isindeq[*vIt] = true;
		}
	}

	if(deq.empty()){
		UG_THROW("topological_ordering_core: Graph is not cycle-free [1]!\n");
	}

	vd_t cur_v;
	size_t k = 0;
	size_t miss = 0;
	while(!deq.empty() && miss < deq.size()){
		cur_v = deq.front();
		deq.pop_front();

		//rotate deque, there is a cycle iff miss==deq.size()
		if(indegs[cur_v] > 0){
			deq.push_back(cur_v);
			++miss;
			continue;
		}

		miss = 0;

		if(inverse){
			o[k++] = cur_v;
		}
		else{
			o[cur_v] = k++;
		}

		for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(cur_v, g); nIt != nEnd; ++nIt){
			--indegs[*nIt];

			if(isindeq[*nIt]){
				continue;
			}

			if(indegs[*nIt] == 0){
				deq.push_front(*nIt);
				isindeq[*nIt] = true;
			}
			else if(indegs[*nIt] > 0){
				deq.push_back(*nIt);
				isindeq[*nIt] = true;
			}
			else{} //ignore vertex
		}
	}

	if(!deq.empty()){
		UG_THROW("topological_ordering_core: Graph is not cycle-free [2]!\n");
	}
}


//for cycle-free matrices only
template <typename TAlgebra, typename O_t>
class TopologicalOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::directedS> G_t;
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

	void compute(){
		topological_ordering_core(o, g, false); //false = no inverse

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

		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){
					boost::add_edge(conn.index(), i, g);
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

	virtual const char* name() const {return "TopologicalOrdering";}

private:
	G_t g;
	O_t o;

	std::deque<vd_t> deq;
	std::vector<int> indegs;
	std::vector<BOOL> isindeq;
};

} //namespace

#endif //guard

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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_OWN_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_OWN_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <set>
#include <algorithm> //reverse
#include <utility> //pair

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"

// is this already implemented somewhere else?
namespace{
template<class T>
double my_abs(T){return 0;}

template<>
double my_abs(double v){return abs(v);}
}

namespace ug{

#ifndef PRINT_GRAPH
#define PRINT_GRAPH
template <typename G_t>
void print_graph(G_t& g){
	typedef typename boost::graph_traits<G_t>::edge_descriptor Edge;

	typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
	for(boost::tie(eIt, eEnd) = boost::edges(g); eIt != eEnd; ++eIt){
		std::pair<Edge, bool> e = boost::edge(boost::source(*eIt, g), boost::target(*eIt, g), g);
		double w = boost::get(boost::edge_weight_t(), g, e.first);
		std::cout << boost::source(*eIt, g) << " -> " << boost::target(*eIt, g) << " ( " << w << " )" << std::endl;
	}
}
#endif


/*

template <typename S_t>
void print(S_t &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << sIt->source << " -> " << sIt->target << " ( " << sIt->w << ")" << std::endl;
	} std::cout << std::endl;
}

template <typename T>
void print(std::set<T> &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << *sIt << " ";
	} std::cout << std::endl;
}

template <typename T>
void print(std::vector<T> &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << *sIt << " ";
	} std::cout << std::endl;
}
*/


template <typename TAlgebra, typename O_t>
class OwnCuthillMcKeeOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::out_edge_iterator oute_iter;

	OwnCuthillMcKeeOrdering() : m_bReverse(false){}

	/// clone constructor
	OwnCuthillMcKeeOrdering( const OwnCuthillMcKeeOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new OwnCuthillMcKeeOrdering<TAlgebra, O_t>(*this));
	}

	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		o.resize(n);

		UG_LOG(name() << "::compute: n = " << n << ", m = " << boost::num_edges(g) << std::endl);
#ifdef DUMP
		std::cout << "graph: " << std::endl; print_graph(g); std::cout << "end graph" << std::endl;
#endif

		std::vector<size_t> indegs(n);

		vIt_t vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			indegs[*vIt] = boost::in_degree(*vIt, g);
		}

		vd s;
		size_t min_s = n;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			if(indegs[*vIt] < min_s){
				s = *vIt;
				min_s = indegs[*vIt];
			}
		}

		UG_LOG("s = " << s << ", deg: " << min_s << "\n");

		for(unsigned i = 0; i < n; ++i){
			o[i] = i;
		}

		g = G_t(0);

		std::cout << "[CuthillMcKeeOrderingWeighted::compute] done. " << std::endl;
	}

	void init(M_t* A, const V_t&){
		init(A);
	}

	void init(M_t* A){
		unsigned n = A->num_rows();

		g = G_t(n);

		for(unsigned i = 0; i < n; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){
					double w;
					w = my_abs(conn.value());
					boost::add_edge(conn.index(), i, w, g);
				}
			}
		}
	}

	void init(M_t* A, const V_t&, const O_t& inv_map){
		init(A, inv_map);
	}

	void init(M_t* A, const O_t& inv_map){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << " on induced matrix of size " << inv_map.size() << "\n");
		#endif
	}

	void check(){
		if(!is_permutation(o)){
			UG_THROW(name() << "::check: Not a permutation!");
		}
	}

	O_t& ordering(){
		return o;
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	virtual const char* name() const {return "OwnCuthillMcKeeOrdering";}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
};

} //namespace


#endif //guard

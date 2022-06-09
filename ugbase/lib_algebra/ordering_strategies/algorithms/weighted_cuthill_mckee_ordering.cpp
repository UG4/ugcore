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

#ifndef __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_Weighted_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_DISC__ORDERING_STRATEGIES_ALGORITHMS_Weighted_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <set>
#include <algorithm> //reverse
#include <utility> //pair
#include <deque>
#include <list>

#include "lib_algebra/ordering_strategies/algorithms/IOrderingAlgorithm.h"
#include "lib_algebra/ordering_strategies/algorithms/util.cpp"

#include <assert.h>
#include "common/error.h"


namespace ug{

namespace{
template<class T>
double my_abs(T v){return 0;}

template<>
double my_abs(double v){
		return abs(v);
}
}


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




#ifndef INDUCED_SUBGRAPH_WEIGHTED
#define INDUCED_SUBGRAPH_WEIGHTED
template <typename G_t, typename M_t>
void induced_subgraph_weighted(G_t& ind_g, M_t* A, const std::vector<size_t>& inv_map){
	size_t n = A->num_rows();
	size_t k = inv_map.size();
	ind_g = G_t(k);

	std::vector<int> ind_map(n, -1);
	for(unsigned i = 0; i < k; ++i){
		ind_map[inv_map[i]] = i;
	}

	typename boost::graph_traits<G_t>::adjacency_iterator nIt, nEnd;
	for(unsigned i = 0; i < inv_map.size(); ++i){
		for(typename M_t::row_iterator conn = A->begin_row(inv_map[i]); conn != A->end_row(inv_map[i]); ++conn){
			if(conn.value() != 0.0 && conn.index() != i){
				int idx = ind_map[conn.index()];
				if(idx >= 0){
					double a = my_abs(conn.value());
					boost::add_edge(idx, i, a, ind_g);
				}
			}
		}
	}
}
#endif

template <typename TAlgebra, typename O_t>
class WeightedCuthillMcKeeOrdering : public IOrderingAlgorithm<TAlgebra, O_t>
{
public:
	typedef typename TAlgebra::matrix_type M_t;
	typedef typename TAlgebra::vector_type V_t;
	typedef IOrderingAlgorithm<TAlgebra, O_t> baseclass;

	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::edge_descriptor ed_t;
	typedef typename boost::graph_traits<G_t>::vertex_iterator vIt_t;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;
	typedef typename boost::graph_traits<G_t>::out_edge_iterator oute_iter;

	WeightedCuthillMcKeeOrdering() : m_bReverse(false){}

	/// clone constructor
	WeightedCuthillMcKeeOrdering( const WeightedCuthillMcKeeOrdering<TAlgebra, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse){}

	SmartPtr<IOrderingAlgorithm<TAlgebra, O_t> > clone()
	{
		return make_sp(new WeightedCuthillMcKeeOrdering<TAlgebra, O_t>(*this));
	}

	inline void unregister_indegree(size_t v, std::vector<size_t>& indegs){
		std::pair<oute_iter, oute_iter> e;
		e = boost::out_edges(v, g);

		for(; e.first != e.second; ++e.first){
			ed_t const& edg = *e.first;
			auto t = boost::target(edg, g);
			--indegs[t];
		 }
	}

	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			UG_THROW(name() << "::compute: Graph is empty!");
			return;
		}

		o.resize(n);

		size_t k = 0;
		std::vector<int> numbering(n, -1);
		std::vector<size_t> indegs(n);
		//std::vector<size_t> front(n);
		std::list<size_t> front;

		vIt_t vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			indegs[*vIt] = boost::in_degree(*vIt, g);
			//source vertex
			if(indegs[*vIt] == 0){
				numbering[*vIt] = k;
				front.push_back(*vIt);
				unregister_indegree(*vIt, indegs);
				++k;
			}
		}
		UG_COND_THROW(k == 0, name() << ": no sources numbered, front empty != n\n");

		UG_LOG("#sources numbered: " << k << "\n");

		std::pair<oute_iter, oute_iter> e;

		//main loop
		for(; k < n;){
			size_t min_indeg_v;
			size_t min_indeg_val = n;

			std::vector<std::list<size_t>::iterator> to_delete;

			for(auto it = front.begin(); it != front.end(); ++it){
				auto candidate = *it;

				size_t misses = 0;
				e = boost::out_edges(candidate, g);
				for(; e.first != e.second; ++e.first){
					ed_t const& edg = *e.first;
					auto t = boost::target(edg, g);

					//unnumbered vertex
					if(numbering[t] < 0){
						if(indegs[t] < min_indeg_val){
							min_indeg_val = indegs[t];
							min_indeg_v = t;
						}
					}
					else{
						++misses;
					}
				 }

				if(misses == boost::out_degree(candidate, g)){
					to_delete.push_back(it);
				}
			}

			for(size_t u = 0; u < to_delete.size(); ++u){
				front.erase(to_delete[u]);
			}

			numbering[min_indeg_v] = k;
			front.push_back(min_indeg_v);

			unregister_indegree(min_indeg_v, indegs);
			++k;
		}

		for(unsigned i = 0; i < n; ++i){
			o[i] = numbering[i];
		}

		g = G_t(0);

		UG_LOG(name() << ": done.\n");
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

	void init(M_t* A, const V_t&, const O_t& inv_map, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
	}

	void init(M_t* A, const O_t& inv_map, const O_t&){
		UG_THROW(name() << "::init: induced subgraph version not implemented yet!");
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

	virtual const char* name() const {return "WeightedCuthillMcKeeOrdering";}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
};

} //namespace


#endif //guard

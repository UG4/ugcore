/*
 * Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_WEIGHTED_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_WEIGHTED_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <set>
#include <algorithm> //reverse
#include <utility> //pair

#include "IOrderingAlgorithm.h"

#include "iters.cpp"
#include "util.cpp"

#include <assert.h>
#include "../../../common/code_marker.h" //error()


namespace ug{

template <typename S_t>
void print(S_t &s){
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << sIt->source << " -> " << sIt->target << " ( " << sIt->w << ")" << std::endl;
	} std::cout << std::endl;
}

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

//#define DUMP

/* TODO

	implement the while loop with iterators

*/


// TODO: rework this if it produces good orderings and needs to be fast..


template <typename M_t, typename O_t>
class WeightedCuthillMcKeeOrdering : public IOrderingAlgorithm<M_t, O_t>
{
public:
	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;

	typedef IOrderingAlgorithm<M_t, O_t> baseclass;

	WeightedCuthillMcKeeOrdering() : m_bReverse(false){}

#if 0
	class next_vertex_iterator{
	public:
		next_vertex_iterator(){}

		bool operator==(const next_vertex_iterator& o) const{}

		bool operator!=(const next_vertex_iterator& o) const{
			return !operator==(o);
		}

		void operator++(){}

		unsigned operator*() const{}
	};


	std::pair<next_vertex_iterator, next_vertex_iterator> make_next_vertex_iterator(G_t* g)
	{
		auto a = next_vertex_iterator(0, g);
		auto b = next_vertex_iterator(boost::num_vertices(g), g);

		return std::make_pair(a, b);
	}
#endif

	double compute_inflow(vd v, std::vector<BOOL>& visited, bool ignore_visited=true){
		double w = .0f;
		typename boost::graph_traits<G_t>::in_edge_iterator in_nIt, in_nEnd;
		for(boost::tie(in_nIt, in_nEnd) = boost::in_edges(v, g); in_nIt != in_nEnd; ++in_nIt){
			if(ignore_visited && visited[v]){ continue; }
			w += abs(boost::get(boost::edge_weight_t(), g, *in_nIt)); //TODO: think about this!
		}
		return w;
	}

	//TODO: do not recompute
	vd min_inflow_vertex(std::set<vd>& front, std::vector<BOOL>& visited){
		double min_w = -1u;
		vd min_vertex;

		assert(front.size());

		typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
		for(auto sIt = front.begin(); sIt != front.end(); ++sIt){
			double w = compute_inflow(*sIt, visited);
#ifdef DUMP
			std::cout << *sIt << ": " << w << std::endl;
#endif
			if(w < min_w){
				min_w = w;
				min_vertex = *sIt;
			}
		}
		return min_vertex;
	}


	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			std::cerr << "graph not set! abort." << std::endl;
			return;
		}

		o.resize(n);
		std::cout << "[CuthillMcKeeOrderingWeighted::compute] n = " << n << ", e = " << boost::num_edges(g) << std::endl;
#ifdef DUMP
		std::cout << "graph: " << std::endl; print_graph(g); std::cout << "end graph" << std::endl;
#endif


		//std::pair<next_vertex_iterator, next_vertex_iterator> tuple_next_vertex_iterators = make_next_vertex_iterator(g);
		//auto nextIt = tuple_next_vertex_iterators.first;
		//auto nextEnd = tuple_next_vertex_iterators.second; TODO: fix end

		vd cur;
		std::set<vd> front;
		std::vector<BOOL> visited(n, false);

		for(unsigned k = 0; k < n; ++k){
			/*
				the algorithm can be tuned here: start vertex and next vertex
			*/

			if(front.size() == 0){
				//initial front - TODO: ordering begins with *boost::vertices(g).first !
				//cur = *unvisited_iterator(visited);
				for(unsigned i = 0; i < boost::num_vertices(g); ++i){
					if(!visited[i]){
						cur = i;
					}
				}
			}
			else{
				cur = min_inflow_vertex(front, visited);
			}
#ifdef DUMP
			std::cout << "next vertex: " << cur << std::endl;
#endif


			visited[cur] = true;

			//insert downstream neighbors to front
			typename boost::graph_traits<G_t>::out_edge_iterator out_nIt, out_nEnd;
			for(boost::tie(out_nIt, out_nEnd) = boost::out_edges(cur, g); out_nIt != out_nEnd; ++out_nIt){
				if(!visited[boost::target(*out_nIt, g)]){
					front.insert(boost::target(*out_nIt, g));
				}
			}

			front.erase(cur);

#ifdef DUMP
			std::cout << "front: "; print(front);
#endif

			o[k] = cur;

#ifdef DUMP
			std::cout << "ordering: ";
			for(unsigned i = 0; i < k; ++i){
				std::cout << o[i] << " ";
			}std::cout << std::endl;
#endif
		}

		if(m_bReverse){
			std::reverse(o.begin(), o.end());
		}

		std::cout << "[CuthillMcKeeOrderingWeighted::compute] done. " << std::endl;
	}

	void check(){
		if(!is_permutation(o)){
			std::cerr << "Not a permutation!" << std::endl;
			print(o);
			error();
		}
	}

	O_t& ordering(){
		return o;
	}

	void set_matrix(M_t* A){
		unsigned rows = A->num_rows();

		g = G_t(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename M_t::row_iterator conn = A->begin_row(i); conn != A->end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					double w;
	#ifdef UG_CPU_1
					w = abs(conn.value()); //TODO: think about this
	#endif
	#ifdef UG_CPU_2
					std::cerr << "[WeightedMatrixGraph] CPU > 1 not implemented yet!" << std::endl;
					error();
	#endif
	#ifdef UG_CPU_3
					std::cerr << "[WeightedMatrixGraph] CPU > 1 not implemented yet!" << std::endl;
					error();
	#endif
					boost::add_edge(i, conn.index(), w, g);
				}
			}
		}
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	std::string config_string() const{
		std::stringstream ss;
		ss << "WeightedCuthillMcKeeOrdering";
		return ss.str();
	}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
};


template <typename M_t, typename O_t>
void weighted_Cuthill_McKee_ordering(M_t &m, bool reverse){
	WeightedCuthillMcKeeOrdering<M_t, O_t> algo();
	algo.set_matrix(m);
	algo.set_reverse(reverse);
	algo.compute();
}

} //namespace


#endif //guard

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

#include <assert.h>
#include "../../../common/code_marker.h" //error()


namespace ug{

//#define DUMP

template <typename S_t>
void print(S_t &s){
#ifdef DUMP
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << sIt->source << " -> " << sIt->target << " ( " << sIt->w << ")" << std::endl;
	} std::cout << std::endl;
#endif
}

template <typename G_t>
void print_graph(G_t &g){
#ifdef DUMP
	typedef typename G_t::edge_descriptor Edge;

	typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
	for(boost::tie(eIt, eEnd) = boost::edges(g); eIt != eEnd; ++eIt){
		std::pair<Edge, bool> e = boost::edge(boost::source(*eIt, g), boost::target(*eIt, g), g);
		double w = boost::get(boost::edge_weight_t(), g, e.first);
		std::cout << boost::source(*eIt, g) << " -> " << boost::target(*eIt, g) << " ( " << w << " )" << std::endl;
	}
#endif
}

template <typename T>
void print(std::set<T> &s){
#ifdef DUMP
	for(auto sIt = s.begin(); sIt != s.end(); ++sIt){
		std::cout << *sIt << " ";
	} std::cout << std::endl;
#endif
}


template <typename T>
void dump(T t){
#ifdef DUMP
	std::cout << t;
#endif
}



/* TODO

	implement the while loop with iterators

*/


// TODO: rework this if it produces good orderings and needs to be fast..


template <typename G_t, typename O_t>
class WeightedCuthillMcKeeOrdering : public IOrderingAlgorithm<O_t>
{
public:
	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;

	WeightedCuthillMcKeeOrdering(G_t &g_in, O_t &o_in, bool reverse) : g(g_in), o(o_in), m_bReverse(reverse){}
	~WeightedCuthillMcKeeOrdering(){}

	class next_vertex_iterator{
	public:
		next_vertex_iterator(vd cur, G_t &g) : _cur(cur), _g(g), _init(true){
			unsigned n = boost::num_vertices(g);
			_visited = std::vector<BOOL>(n, false);
			if(_cur < n){ //TODO: should be in make_next_...
				operator++();
			}
			_k = 0;
		}

		bool operator==(const next_vertex_iterator& o) const{
			return _cur==o._cur;
		}

		bool operator!=(const next_vertex_iterator& o) const{
			return !operator==(o);
		}



		double compute_inflow(vd v, bool ignore_visited=true){
			double w = .0f;
			typename boost::graph_traits<G_t>::in_edge_iterator in_nIt, in_nEnd;
			for(boost::tie(in_nIt, in_nEnd) = boost::in_edges(v, _g); in_nIt != in_nEnd; ++in_nIt){
				if(ignore_visited && _visited[v]){ continue; }
				w += abs(boost::get(boost::edge_weight_t(), _g, *in_nIt)); //TODO: think about this!
			}
			return w;
		}


		//TODO: do not recompute
		vd min_inflow_vertex(){
			double min_w = -1u;
			vd min_vertex;

			assert(_front.size());

			typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
			for(auto sIt = _front.begin(); sIt != _front.end(); ++sIt){
				double w = compute_inflow(*sIt);
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

		void operator++(){
			++_k;
			//TODO
			if(_k == boost::num_vertices(_g)){
				return;
			}


			/*
				the algorithm can be tuned here: start vertex and next vertex
			*/

			if(_front.size() == 0){
				//initial front - TODO: ordering begins with *boost::vertices(g).first !
				_cur = *unvisited_iterator(_visited); //TODO: see initial front below
			}
			else{
				_cur = min_inflow_vertex();
			}

#ifdef DUMP
			std::cout << "next vertex: " << _cur << std::endl;
#endif

			_visited[_cur] = true;

			//insert downstream neighbors to front
			typename boost::graph_traits<G_t>::out_edge_iterator out_nIt, out_nEnd;
			for(boost::tie(out_nIt, out_nEnd) = boost::out_edges(_cur, _g); out_nIt != out_nEnd; ++out_nIt){
				if(!_visited[boost::target(*out_nIt, _g)]){
					_front.insert(boost::target(*out_nIt, _g));
				}
			}

			_front.erase(_cur);

#ifdef DUMP
			std::cout << "front: "; print(_front);
#endif
		}

		unsigned operator*() const{
			return _cur;
		}

	private:
		unsigned _k;
		vd _cur;
		std::set<vd> _front;
		std::vector<BOOL> _visited;
		G_t& _g;
		bool _init;
	};


	std::pair<next_vertex_iterator, next_vertex_iterator> make_next_vertex_iterator(G_t &g)
	{
		auto a = next_vertex_iterator(0, g);
		auto b = next_vertex_iterator(boost::num_vertices(g), g);

		return std::make_pair(a, b);
	}


	//overload
	void compute(){
#ifdef DUMP
		std::cout << "[CuthillMcKeeOrderingWeighted::compute]" << std::endl;
#endif

		unsigned n = boost::num_vertices(g);
		o.resize(n);

#ifdef DUMP
		std::cout << "graph: " << std::endl; print_graph(g);;
#endif


		std::pair<next_vertex_iterator, next_vertex_iterator> tuple_next_vertex_iterators = make_next_vertex_iterator(g);
		auto nextIt = tuple_next_vertex_iterators.first;
		//auto nextEnd = tuple_next_vertex_iterators.second; TODO: fix end
		unsigned k = 0;

		for(; k < n; ++nextIt){
			o[k++] = *nextIt;

#ifdef DUMP
			std::cout << "ordering: ";
			for(unsigned i = 0; i < n; ++i){
				std::cout << o[i] << " ";
			}std::cout << std::endl;
#endif
		}

		if(m_bReverse){
			std::reverse(o.begin(), o.end());
		}
	}

	O_t& ordering(){
		return o;
	}

private:
	G_t& g;
	O_t& o;

	bool m_bReverse;
};


template <typename G_t, typename O_t>
void weighted_Cuthill_McKee_ordering(G_t &g, O_t &o, bool reverse){
	WeightedCuthillMcKeeOrdering<G_t, O_t> algo(g, o, reverse);
	algo.compute();
}

} //namespace


#endif //guard

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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_WITTUM_DOWNWIND_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_WITTUM_DOWNWIND_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <utility> //pair

#include "IOrderingAlgorithm.h"

#include "util.cpp"

#include <assert.h>
#include "../../../common/code_marker.h" //error()

#include <boost/graph/strong_components.hpp>


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

//for sorting
struct Bla{
	size_t v;
	unsigned f;
};

bool compBla(Bla a, Bla b){
	return a.f > b.f;
}


class StronglyConnectedComponents
{
public:
	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;

	StronglyConnectedComponents(G_t &g) : _g(g){}

	void create_reverse(){
		_g_rev = G_t(boost::num_vertices(_g));

		typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
		for(boost::tie(eIt, eEnd) = boost::edges(_g); eIt != eEnd; ++eIt){
			double w = abs(boost::get(boost::edge_weight_t(), _g, *eIt));

			boost::add_edge(boost::target(*eIt, _g), boost::source(*eIt, _g), w, _g_rev);
		}
	}

	void DFS_loop_subroutine1(){
		_t = 0;
		_s = -1;

		for(int i = boost::num_vertices(_g_rev)-1; i >= 0; --i){
			if(!visited[i]){
				_s = i;
				DFS(i, _g_rev);
			}
		}
	}

	void DFS_loop_subroutine2(){
		_t = 0;
		_s = -1;

		for(unsigned i = 0; i < boost::num_vertices(_g); ++i){
			if(!visited[decreasing[i]]){
				_s = i;
				DFS(decreasing[i], _g);
			}
		}
	}

	void DFS(vd v, G_t &g){
		visited[v] = true;
		leader[v] = _s;

		adj_iter nIt, nEnd;
		for(boost::tie(nIt, nEnd) = boost::adjacent_vertices(v, g); nIt != nEnd; ++nIt){
			if(!visited[*nIt]){
				DFS(*nIt, g);
			}

			++_t;
			finished[v] = _t;
		}
	}

	void do_it(){
		create_reverse();

		size_t n = boost::num_vertices(_g);
		visited = std::vector<BOOL>(n, false);
		leader = std::vector<vd>(n, -1);
		finished = std::vector<unsigned>(n, -1);
		decreasing = std::vector<vd>(n, -1);

		std::vector<Bla> bla(n);
		for(unsigned i = 0; i < n; ++i){
			bla[i].v = i;
			bla[i].f = finished[i];
		}

		//sort 'decreasing' according to 'finished'
		std::sort(bla.begin(), bla.end(), compBla);

		for(unsigned i = 0; i < n; ++i){
			decreasing[i] = bla[i].v;
		}

		DFS_loop_subroutine1();
		visited = std::vector<BOOL>(n, false);
		DFS_loop_subroutine2();
	}
private:
	G_t &_g;
	G_t _g_rev;

	unsigned _t;
	vd _s;

	std::vector<BOOL> visited;
	std::vector<vd> leader;
	std::vector<unsigned> finished;
	std::vector<vd> decreasing;
};


template <typename M_t, typename O_t>
class WittumDownwindOrdering : public IOrderingAlgorithm<M_t, O_t>
{
public:
	typedef boost::property<boost::edge_weight_t, double> EdgeWeightProperty;
	typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS, boost::no_property, EdgeWeightProperty> G_t;

	typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
	typedef typename boost::graph_traits<G_t>::adjacency_iterator adj_iter;

	typedef IOrderingAlgorithm<M_t, O_t> baseclass;

	WittumDownwindOrdering(){}

	/// clone constructor
	WittumDownwindOrdering( const WittumDownwindOrdering<M_t, O_t> &parent )
			: baseclass(){}

	SmartPtr<IOrderingAlgorithm<M_t, O_t> > clone()
	{
		return make_sp(new WittumDownwindOrdering<M_t, O_t>(*this));
	}

	size_t compute_inedges(vd v){
		size_t k = 0;
		typename boost::graph_traits<G_t>::in_edge_iterator in_nIt, in_nEnd;
		for(boost::tie(in_nIt, in_nEnd) = boost::in_edges(v, g); in_nIt != in_nEnd; ++in_nIt){
			if(!visited[boost::source(*in_nIt, g)]){
				++k;
			}
		}

		return k;
	}

	double compute_inflow(vd v){
		double w = .0f;
		typename boost::graph_traits<G_t>::in_edge_iterator in_nIt, in_nEnd;
		for(boost::tie(in_nIt, in_nEnd) = boost::in_edges(v, g); in_nIt != in_nEnd; ++in_nIt){
			w += abs(boost::get(boost::edge_weight_t(), g, *in_nIt));
		}
		return w;
	}

	//TODO: do not recompute
	vd min_inedges_vertex(){
		double min_k = -1u;
		vd min_vertex;

		typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; vIt++){
			if(visited[*vIt]){ continue; }
			size_t indeg = compute_inedges(*vIt);
			if(indeg < min_k){
				min_k = indeg;
				min_vertex = *vIt;
			}
		}

		return min_vertex;
	}

	void remove_smallest_edge(){
		double min_w = 10e9;
		typename boost::graph_traits<G_t>::edge_descriptor min_e;
		typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
		for(boost::tie(eIt, eEnd) = boost::edges(g); eIt != eEnd; ++eIt){
			double w = abs(boost::get(boost::edge_weight_t(), g, *eIt));

			if(w < min_w){
				min_w = w;
				min_e = *eIt;
			}
		}

		std::cout << "remove edge with weight " << min_w << ", e = " << boost::num_edges(g) << std::endl;
		boost::remove_edge(min_e, g);
	}

	void strongly_connected_components(){
		std::cout << "strongly_connected_components()" << std::endl;
		size_t n = boost::num_vertices(g);
 		std::vector<int> component(n);
		std::vector<int> discover_time(n);
    		std::vector<boost::default_color_type> color(n);
    		std::vector<vd> root(n);
    		int num = boost::strong_components(g,
			boost::make_iterator_property_map(component.begin(), boost::get(boost::vertex_index, g)),
			boost::root_map(boost::make_iterator_property_map(root.begin(), boost::get(boost::vertex_index, g)))
			    .color_map(
				boost::make_iterator_property_map(color.begin(), boost::get(boost::vertex_index, g)))
			    .discover_time_map(boost::make_iterator_property_map(
				discover_time.begin(), boost::get(boost::vertex_index, g))));

		std::vector<size_t> comp_sizes(num, 0);

		for(unsigned i = 0; i < n; ++i){
			++comp_sizes[component[i]];
		}

		std::cout << num << " components" << std::endl;

		for(unsigned i = 0; i < num; ++i){
			std::cout << "comp " << i << ", # " << comp_sizes[i] << std::endl;
		}

		for(unsigned i = 0; i < n; ++i){
			std::cout << i << ": " << component[i] << std::endl;
		}
	}

	//overload
	void compute(){
		unsigned n = boost::num_vertices(g);

		if(n == 0){
			std::cerr << "graph not set! abort." << std::endl;
			return;
		}

		visited = std::vector<BOOL>(n, false);
		o.resize(n);

		std::cout << "[WittumDownwindOrdering::compute] n = " << n << ", e = " << boost::num_edges(g) << std::endl;
		print_graph(g);

		for(unsigned i = 0; i < n; ++i){
			vd v = min_inedges_vertex();

			size_t deg = compute_inedges(v);

			while(deg > 0){
				strongly_connected_components();
				remove_smallest_edge();
				v = min_inedges_vertex();
				deg = compute_inedges(v);
			}

			std::cout << "i=" << i << ", indegree " << deg << std::endl;
			o[i] = v;
			visited[v] = true;
		}

		std::cout << "ordering: ";
		for(unsigned i = 0; i < n; ++i){
			std::cout << o[i] << " ";
		} std::cout << std::endl;

		g = G_t(0);

		std::cout << "[WittumDownwindOrdering::compute] done. " << std::endl;

		check();
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

	std::string config_string() const{
		std::stringstream ss;
		ss << "WittumDownwindOrdering";
		return ss.str();
	}

private:
	G_t g;
	O_t o;

	std::vector<BOOL> visited;
};

} //namespace


#endif //guard

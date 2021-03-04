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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_SHORTEST_PATHS_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_SHORTEST_PATHS_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/dijkstra_shortest_paths.hpp>

#include "IOrderingAlgorithm.h"

#include "../execution/util.cpp"

namespace ug{

//for sorting
struct Blo{
	size_t v;
	double w;
};

bool compBlo(Blo a, Blo b){
	return a.w < b.w;
}


template <typename M_t, typename G_t, typename O_t>
class BoostShortestPathsOrdering : public IOrderingAlgorithm<M_t, G_t, O_t>
{
public:
	typedef IOrderingAlgorithm<M_t, G_t, O_t> baseclass;
	typedef typename baseclass::Type Type;

	//BoostShortestPathsOrdering(G_t &g_in, O_t &o_in) : g(&g_in), o(&o_in), own_o(false){}
	BoostShortestPathsOrdering() : own_o(false){}
	~BoostShortestPathsOrdering(){
		if(own_o){ delete o; }
	}

	void compute(){
		if(!g){
			std::cerr << "graph not set! abort." << std::endl;
			return;
		}

		if(!o){
			own_o = true;
			o = new O_t;
		}

		typedef typename boost::graph_traits<G_t>::vertex_descriptor vd;
		std::vector<vd> p(boost::num_vertices(*g)); //parents
		std::vector<int> d(num_vertices(*g)); //distances

		vd s = *boost::vertices(*g).first; //start vertex	//TODO: choose a vertex strategically

		boost::dijkstra_shortest_paths(*g, s, boost::predecessor_map(&p[0]).distance_map(&d[0]));

		std::vector<Blo> blo(boost::num_vertices(*g));
		for(unsigned i = 0; i < boost::num_vertices(*g); ++i){
			blo[i].v = i;
			blo[i].w = d[i];
		}

		//sort o according to d
		std::sort(blo.begin(), blo.end(), compBlo);

		o->resize(boost::num_vertices(*g));
		for(unsigned i = 0; i < boost::num_vertices(*g); ++i){
			(*o)[i] = blo[i].v;
		}
	}

	void check(){
		if(!is_permutation(*o)){
			std::cerr << "Not a permutation!" << std::endl;
			error();
		}
	}

	O_t* ordering(){
		return o;
	}

	const Type type(){
		return mytype;
	} 

	void set_graph(G_t* graph){
		g = graph;
	}

	void set_matrix(M_t*){}

	void set_ordering(O_t* ordering){
		o = ordering;
	}
private:
	G_t* g;
	O_t* o;

	bool own_o;

	static const Type mytype = Type::GRAPH_BASED;
};


template <typename M_t, typename G_t, typename O_t>
void boost_shortest_paths_ordering(G_t &g, O_t &o){
	BoostShortestPathsOrdering<M_t, G_t, O_t> algo(&g, &o);
	algo.compute();
}


} //namespace

#endif //guard

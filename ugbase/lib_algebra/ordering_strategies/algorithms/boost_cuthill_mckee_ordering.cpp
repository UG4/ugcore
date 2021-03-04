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
 
#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_ORDERING__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_BOOST_CUTHILL_MCKEE_ORDERING__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include <boost/graph/cuthill_mckee_ordering.hpp>

#include "IOrderingAlgorithm.h"

#include "../execution/util.cpp"

namespace ug{


/* boost graph type for Cuthill-McKee */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t,
			 boost::default_color_type,
			 boost::property<boost::vertex_degree_t, int> > >
				Graph_t;


template <typename G_t>
void copy_graph_for_boost_Cuthill_McKee(G_t &orig, Graph_t &copy){
	copy = Graph_t(boost::num_vertices(orig));

	typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
	for(boost::tie(eIt, eEnd) = boost::edges(orig); eIt != eEnd; ++eIt){
			boost::add_edge(boost::source(*eIt, orig), boost::target(*eIt, orig), copy);
	}
}


template <typename M_t, typename G_t, typename O_t>
class BoostCuthillMcKeeOrdering : public IOrderingAlgorithm<M_t, G_t, O_t>
{
public:
	typedef IOrderingAlgorithm<M_t, G_t, O_t> baseclass;
	typedef typename baseclass::Type Type;
	BoostCuthillMcKeeOrdering() : m_bReverse(false), own_o(false){}
	~BoostCuthillMcKeeOrdering(){
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

		typedef boost::graph_traits<Graph_t>::vertex_descriptor Vertex_t;
		typedef boost::graph_traits<Graph_t>::vertices_size_type size_type;

		Graph_t h;
		copy_graph_for_boost_Cuthill_McKee(*g, h); //explicit copy

		boost::property_map<Graph_t, boost::vertex_degree_t>::type deg = get(boost::vertex_degree, h);
		boost::graph_traits<Graph_t>::vertex_iterator vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(h); vIt != vEnd; ++vIt){
			deg[*vIt] = boost::degree(*vIt, h);
		}

		boost::property_map<Graph_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, h);

		std::vector<Vertex_t> inv_perm(boost::num_vertices(h));

		if(m_bReverse){
			boost::cuthill_mckee_ordering(h, inv_perm.rbegin(), get(boost::vertex_color, h), boost::make_degree_map(h));
		}
		else{
			boost::cuthill_mckee_ordering(h, inv_perm.begin(), get(boost::vertex_color, h), boost::make_degree_map(h));
		}

		o->resize(boost::num_vertices(*g));

		for(unsigned i = 0; i != inv_perm.size(); ++i){
			(*o)[index_map[inv_perm[i]]] = i;
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

	void set_reverse(bool b){
		m_bReverse = b;
	}

private:
	G_t* g;
	O_t* o;

	bool m_bReverse;

	bool own_o;

	static const Type mytype = Type::GRAPH_BASED;
};


template <typename M_t, typename G_t, typename O_t>
void boost_Cuthill_McKee_ordering(G_t &g, O_t &o, bool reverse){
	BoostCuthillMcKeeOrdering<M_t, G_t, O_t> algo(g, o, reverse);
	algo.compute();
}


} //namespace

#endif //guard

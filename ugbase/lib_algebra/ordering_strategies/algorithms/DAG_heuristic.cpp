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

#ifndef __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_DAG_HEURISTIC__
#define __UG__LIB_ALGEBRA__ORDERING_STRATEGIES_ALGORITHMS_DAG_HEURISTIC__

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

//debug
#include "common/error.h"
#include "common/log.h"

namespace ug{

template <typename G_t, typename O_t>
class DAGHeuristic
{
public:
	//typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::bidirectionalS> G_t;

	DAGHeuristic(){}

	/// clone constructor
	DAGHeuristic( const DAGHeuristic<G_t, O_t> &parent ){}

	SmartPtr<DAGHeuristic<G_t, O_t> > clone()
	{
		return make_sp(new DAGHeuristic<G_t, O_t>(*this));
	}

	void compute(){
		typename boost::graph_traits<G_t>::vertex_descriptor s, t;
		typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
		for(boost::tie(eIt, eEnd) = boost::edges(*g); eIt != eEnd; ++eIt){
			s = boost::source(*eIt, *g);
			t = boost::target(*eIt, *g);
			if((*o)[s] < (*o)[t]){
				boost::add_edge((*o)[s], (*o)[t], g1);
			}
			else{
				boost::add_edge((*o)[s], (*o)[t], g2);
			}
		}

		std::cout << "edges g1: " << boost::num_edges(g1) << std::endl;
		std::cout << "edges g2: " << boost::num_edges(g2) << std::endl;
		std::cout << "edges total: " << boost::num_edges(*g) << std::endl;

		g1 = G_t(0);
		g2 = G_t(0);

		#ifdef UG_DEBUG
		check();
		#endif
	}

	void init(G_t& g_in, O_t& o_in){
		//TODO: replace this by UG_DLOG if permutation_util does not depend on this file anymore
		#ifdef UG_ENABLE_DEBUG_LOGS
		UG_LOG("Using " << name() << "\n");
		#endif

		g = &g_in;
		o = &o_in;

		size_t n = boost::num_vertices(g_in);
		g1 = G_t(n);
		g2 = G_t(n);
	}

	virtual const char* name() const {return "DAGHeuristic";}

private:
	G_t* g;
	G_t g1, g2;
	O_t* o;
};

} //namespace

#endif //guard

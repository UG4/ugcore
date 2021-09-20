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
#include "util.cpp"

#include "common/code_marker.h"

namespace ug{


/* boost graph type for Cuthill-McKee */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS,
	boost::property<boost::vertex_color_t,
			 boost::default_color_type,
			 boost::property<boost::vertex_degree_t, int> > >
				Graph_t;


template <typename G_t>
void print_graph_unweighted(G_t& g){
	typedef typename boost::graph_traits<G_t>::edge_descriptor Edge;

	typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
	for(boost::tie(eIt, eEnd) = boost::edges(g); eIt != eEnd; ++eIt){
		std::cout << boost::source(*eIt, g) << " -> " << boost::target(*eIt, g) << std::endl;
	}
}


template <typename M_t, typename O_t=std::vector<size_t> >
class BoostCuthillMcKeeOrdering : public IOrderingAlgorithm<M_t, O_t>
{
public:
	typedef Graph_t G_t;
	typedef IOrderingAlgorithm<M_t, O_t> baseclass;

	BoostCuthillMcKeeOrdering() : m_bReverse(false){}

	/// clone constructor
	BoostCuthillMcKeeOrdering( const BoostCuthillMcKeeOrdering<M_t, O_t> &parent )
			: baseclass(), m_bReverse(parent.m_bReverse){}

	SmartPtr<IOrderingAlgorithm<M_t, O_t> > clone()
	{
		return make_sp(new BoostCuthillMcKeeOrdering<M_t, O_t>(*this));
	}

	void compute(){
		//std::cout << "graph: " << std::endl; print_graph_unweighted(g); std::cout << "end graph" << std::endl;

		unsigned n = boost::num_vertices(g);

		if(n == 0){
			std::cerr << "graph not set! abort." << std::endl;
			return;
		}

		boost::property_map<Graph_t, boost::vertex_degree_t>::type deg = get(boost::vertex_degree, g);
		boost::graph_traits<Graph_t>::vertex_iterator vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(g); vIt != vEnd; ++vIt){
			deg[*vIt] = boost::degree(*vIt, g);
		}

		boost::property_map<Graph_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, g);

		typedef boost::graph_traits<Graph_t>::vertex_descriptor Vertex_t;
		std::vector<Vertex_t> inv_perm(boost::num_vertices(g));

		if(m_bReverse){
			boost::cuthill_mckee_ordering(g, inv_perm.rbegin(), get(boost::vertex_color, g), boost::make_degree_map(g));
		}
		else{
			boost::cuthill_mckee_ordering(g, inv_perm.begin(), get(boost::vertex_color, g), boost::make_degree_map(g));
		}

		o.resize(boost::num_vertices(g));

		for(unsigned i = 0; i != inv_perm.size(); ++i){
			o[index_map[inv_perm[i]]] = i;
		}

		g = G_t(0);
	}

	void check(){
		if(!is_permutation(o)){
			std::cerr << "Not a permutation!" << std::endl;
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
					boost::add_edge(i, conn.index(), g);
				}
			}
		}
	}

	void set_reverse(bool b){
		m_bReverse = b;
	}

	std::string config_string() const{
		std::stringstream ss;
		ss << "BoostCuthillMcKeeOrdering";
		return ss.str();
	}

private:
	G_t g;
	O_t o;

	bool m_bReverse;
};


template <typename M_t, typename O_t>
void boost_Cuthill_McKee_ordering(M_t &m, bool reverse){
	BoostCuthillMcKeeOrdering<M_t, O_t> algo();
	algo.set_reverse(reverse);
	algo.set_matrix(m);
	algo.compute();
}


} //namespace

#endif //guard

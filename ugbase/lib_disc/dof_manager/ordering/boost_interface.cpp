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
 
#ifndef __H__UG__LIB_DISC__DOF_MANAGER__BOOST_INTERFACE__CPP
#define __H__UG__LIB_DISC__DOF_MANAGER__BOOST_INTERFACE__CPP

#include <chrono>

#include "lib_disc/domain.h"
#include "lib_disc/function_spaces/approximation_space.h"

#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/graph_traits.hpp>
#include <boost/graph/properties.hpp>

#include "boost_interface.hpp"

#include "types.hpp"
#include "util.cpp"
#include "graph_assembly.cpp"

#include <boost/graph/cuthill_mckee_ordering.hpp>
#include <boost/graph/minimum_degree_ordering.hpp>

namespace ug{


/* boost graph type for Cuthill-McKee */
typedef boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, 
	boost::property<boost::vertex_color_t, 
			 boost::default_color_type, 
			 boost::property<boost::vertex_degree_t, int> > > 
				Graph_t;



/* Cuthill-McKee */

template <typename G_t, typename O_t>
void boost_Cuthill_McKee_ordering(G_t &g, O_t &o, bool reverse){
	typedef boost::graph_traits<Graph_t>::vertex_descriptor Vertex_t;
	typedef boost::graph_traits<Graph_t>::vertices_size_type size_type;
	
	Graph_t h;
	copy_graph(g, h);

	boost::property_map<Graph_t, boost::vertex_degree_t>::type deg = get(boost::vertex_degree, h);
	boost::graph_traits<Graph_t>::vertex_iterator vIt, vEnd;
	for(boost::tie(vIt, vEnd) = boost::vertices(h); vIt != vEnd; ++vIt){
    		deg[*vIt] = boost::degree(*vIt, h);
    	}

  	boost::property_map<Graph_t, boost::vertex_index_t>::type index_map = get(boost::vertex_index, h);

	std::vector<Vertex_t> inv_perm(boost::num_vertices(h));

        if(reverse){
        	boost::cuthill_mckee_ordering(h, inv_perm.rbegin(), get(boost::vertex_color, h), boost::make_degree_map(h));
        }
        else{
        	boost::cuthill_mckee_ordering(h, inv_perm.begin(), get(boost::vertex_color, h), boost::make_degree_map(h));
        }
                           
	o.resize(boost::num_vertices(g));
	
	for(unsigned i = 0; i != inv_perm.size(); ++i){
		o[index_map[inv_perm[i]]] = i;
	}
}



/* minimum degree */

  //Important Note: This implementation requires the BGL graph to be
  //directed.  Therefore, nonzero entry (i, j) in a symmetrical matrix
  //A coresponds to two directed edges (i->j and j->i).
template <typename G_t, typename O_t>
void boost_minimum_degree_ordering(G_t &G, O_t &O, O_t &iO){
    unsigned n = boost::num_vertices(G);
    unsigned e = boost::num_edges(G);

    O.resize(n);
    unsigned i = 0;

    if(n == 0){
    	error();
    }
    
    else if(n*(n-1u)==e || e==0){	error();
        typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
        for(boost::tie(vIt, vEnd) = boost::vertices(G); vIt != vEnd; vIt++){
            O[i++] = *vIt;
        }
    }

    std::vector<int> inverse_perm(n, 0);
    std::vector<int> supernode_sizes(n, 1);
    auto id = boost::get(boost::vertex_index, G);
    std::vector<int> degree(n, 0);

    /*
     * (Graph& G,
     *  DegreeMap degree,
     *  InversePermutationMap inverse_perm,
     *  PermutationMap perm,
     *  SuperNodeMap supernode_size,
     *  int delta,
     *  VertexIndexMap vertex_index_map)
     */

    boost::minimum_degree_ordering
             (G,
              boost::make_iterator_property_map(&degree[0], id, degree[0]),
              &iO[0],
              &O[0],
              boost::make_iterator_property_map(&supernode_sizes[0], id, supernode_sizes[0]),
              0,
              id
              );
}

template <typename G_t, typename O_t>
void boost_minimum_degree_ordering(G_t &G, O_t &O){
    O_t iO(boost::num_vertices(G), 0);
    return boost_minimum_degree_ordering(G, O, iO);
}

// -----------------------------


template <typename TDomain>
void OrderBoostCuthillMcKee(SmartPtr<ApproximationSpace<TDomain> > approxSpace, bool reverse){
	unsigned _n = approxSpace->dof_distributions().size();
	
	auto start = std::chrono::high_resolution_clock::now();

	for(unsigned i = 0; i < _n; ++i){
		graph_bidirectional_t g;
		make_boost_graph_ug4(g, approxSpace, i);
		
		std::vector<size_t> o;
		boost_Cuthill_McKee_ordering(g, o, reverse);	

		if(!is_permutation(o)){
			std::cout << "[OrderBoostCuthillMcKee] Not a permutation!" << std::endl;
			error();
		}
		
		approxSpace->dof_distributions()[i]->permute_indices(o);
    	}
    	
    	auto end = std::chrono::high_resolution_clock::now();
    	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    	std::cout << "[OrderBoostCuthillMcKee] took " << elapsed.count() << " ms." << std::endl;
}

template <typename TDomain>
void OrderBoostMinimumDegree(SmartPtr<ApproximationSpace<TDomain> > approxSpace){
	unsigned _n = approxSpace->dof_distributions().size();

	auto start = std::chrono::high_resolution_clock::now();

	for(unsigned i = 0; i < _n; ++i){
		graph_bidirectional_t g;
		make_boost_graph_ug4(g, approxSpace, i);
		
		std::vector<size_t> o;
		boost_minimum_degree_ordering(g, o);	

		if(!is_permutation(o)){
			std::cout << "[OrderBoostMinimumDegree] Not a permutation!" << std::endl;
			error();
		}
		
		approxSpace->dof_distributions()[i]->permute_indices(o);	
    	}

	auto end = std::chrono::high_resolution_clock::now();
    	auto elapsed = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
    	std::cout << "[OrderBoostMinimumDegree] took " << elapsed.count() << " ms." << std::endl;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__DOF_MANAGER__BOOST_INTERFACE__CPP */



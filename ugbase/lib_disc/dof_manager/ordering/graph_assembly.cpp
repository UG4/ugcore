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
 
#ifndef __H__UG__LIB_DISC__DOF_MANAGER__GRAPH_ASSEMBLY__HPP
#define __H__UG__LIB_DISC__DOF_MANAGER__GRAPH_ASSEMBLY__HPP

#include "lib_disc/function_spaces/approximation_space.h"


/* create graph from grid. */

template <typename TDomain, typename G_t>
void make_boost_graph_grid(G_t &G, SmartPtr<ug::ApproximationSpace<TDomain> > approxSpace, size_t grid_level){
	std::vector<SmartPtr<ug::DoFDistribution> > vDD = approxSpace->dof_distributions();
	auto dofDistr = vDD[grid_level];

	std::vector<std::vector<size_t> > vvConnection;

	unsigned num_vert = dofDistr->num_indices();
	vvConnection.resize(num_vert);
	try{
		dofDistr->get_connections<ug::Edge>(vvConnection); //this is the original grid
	}
	UG_CATCH_THROW("[make_boost_graph_grid] No adjacency graph available.");

	G = G_t(vvConnection.size());

	for(unsigned i = 0; i < vvConnection.size(); ++i){
		for(unsigned j = 1; j < vvConnection[i].size(); ++j){ //first element is the index of the vertex itself
			if(i != vvConnection[i][j]){ //no self loops!
				boost::add_edge(i, vvConnection[i][j], G);
				//boost::add_edge(vvConnection[i][j], i, G);
			}
		}
	}
}

template <typename TDomain, typename G_t>
void make_boost_graph_ug4(G_t &G, SmartPtr<ug::ApproximationSpace<TDomain> > approxSpace, size_t grid_level){
	std::vector<SmartPtr<ug::DoFDistribution> > vDD = approxSpace->dof_distributions();
	auto dofDistr = vDD[grid_level];

	std::vector<std::vector<size_t> > vvConnection;
	try{
		dofDistr->get_connections(vvConnection); //this the element-graph
	}
	UG_CATCH_THROW("[make_boost_graph_ug4] No adjacency graph available.");

	G = G_t(vvConnection.size());

	for(unsigned i = 0; i < vvConnection.size(); ++i){
		for(unsigned j = 1; j < vvConnection[i].size(); ++j){ //first element is the index of the vertex itself
			if(i != vvConnection[i][j]){ //no self loops!
				boost::add_edge(i, vvConnection[i][j], G);
				//boost::add_edge(vvConnection[i][j], i, G);
			}
		}
	}
}


/* create graph from matrix. */

template <typename TAlgebra, typename G_t>
void make_boost_graph_from_matrix_ug4(G_t &G, typename TAlgebra::matrix_type &A){
	std::cout << "[make_boost_graph_from_matrix_ug4]" << std::endl;
	unsigned rows = A.num_rows();

	G = G_t(A.num_rows());

	std::cout << "added vertices: " << boost::num_vertices(G) << std::endl;


	for(unsigned i = 0; i < rows; i++){
		for(typename TAlgebra::matrix_type::row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn){
    			std::cout << conn.value() << " ";
			if(conn.value() != 0.0){
				boost::add_edge(i, conn.index(), G);
			}
		} std::cout << std::endl;
	}

	std::cout << "[make_boost_graph_from_matrix_ug4] " <<  "graph has " << boost::num_vertices(G) << " vertices and " << boost::num_edges(G) << " edges" << std::endl;
}


/* create weighted graph from matrix. */

template <typename TAlgebra, typename G_t>
void make_boost_graph_from_matrix_weighted_ug4(G_t &G, typename TAlgebra::matrix_type &A){
	std::cout << "[make_boost_graph_from_matrix_weighted_ug4]" << std::endl;
	unsigned rows = A.num_rows();
	unsigned cols = A.num_cols();

	G = G_t(A.num_rows());

	std::cout << "added vertices: " << boost::num_vertices(G) << std::endl;


	typedef typename TAlgebra::matrix_type matrix_type;
	typedef typename matrix_type::value_type block_type;

	typedef typename TAlgebra::vector_type vector_type;
	typedef typename vector_type::value_type smallvec_type;


	//for(typename TMatrixType::const_row_iterator connij = A.begin_row(i); connij != A.end_row(i); ++connij)

	for(unsigned i = 0; i < rows; i++){
		for(unsigned j = 0; j < cols; j++){
    			block_type b = A(i, j);

			std::cout << b(0,0) << " ";

			if(b != 0.0){
#ifdef UG_CPU_1
				double w = b;
#endif
#ifdef UG_CPU_2
				//double w = b;
#endif
#ifdef UG_CPU_3
				//double w = b;
#endif
				//boost::add_edge(i, conn.index(), w, G);
			}
		} std::cout << std::endl;
	}

	std::cout << "[make_boost_graph_from_matrix_weighted_ug4] " <<  "graph has " << boost::num_vertices(G) << " vertices and " << boost::num_edges(G) << " edges" << std::endl;
}

#endif

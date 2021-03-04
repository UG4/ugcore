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


#ifndef __ORDERING_STRATEGIES_METAGRAPH_WEIGHTED_CLIQUE_MATRIX_GRAPH__
#define __ORDERING_STRATEGIES_METAGRAPH_WEIGHTED_CLIQUE_MATRIX_GRAPH__

#include <boost/graph/graph_traits.hpp>

#include "IMetaGraph.h"

#include "../../../common/code_marker.h" //error()

namespace ug{


/*
	Creates a graph g of type G_t based on a matrix A and accessible via graph().
	The graph g is defines as follows:
		V(g) = {0, .., numrow(A)-1}
		E(g) = {{i, j} : A(i, j) != 0 or A(i, k) != 0 and A(k, j) != 0 for some k}
		w({i, j}) = A(i, j)   			| A(i, j) != 0
			     A(i, k) + A(k, j)         | A(i, k) !=0 and A(k, j) !+ 0 for some k. 
*/

/* TODO: currently for #ifdef UG_CPU_1# only. generalize this. */
template <typename G_t, typename matrix_type>
class WeightedCliqueMatrixGraph : public IMetaGraph<G_t>{
public:
	typedef IMetaGraph<G_t> baseclass;

	WeightedCliqueMatrixGraph(matrix_type& A){
		unsigned rows = A.num_rows();

		G_t h(rows);

		for(unsigned i = 0; i < rows; i++){
			for(typename matrix_type::row_iterator conn = A.begin_row(i); conn != A.end_row(i); ++conn){
				if(conn.value() != 0.0 && conn.index() != i){ //TODO: think about this!!
					double w;
	#ifdef UG_CPU_1
					w = abs(conn.value()); //TODO: think about this
	#endif
	#ifdef UG_CPU_2
					std::cerr << "[WeightedCliqueMatrixGraph] CPU > 1 not implemented yet!" << std::endl;
					error();
	#endif
	#ifdef UG_CPU_3
					std::cerr << "[WeightedCliqueMatrixGraph] CPU > 1 not implemented yet!" << std::endl;
					error();
	#endif

					if(w != 0){
						boost::add_edge(i, conn.index(), w, h);
					}
				}
			}
		}

		baseclass::g = G_t(rows);

		typedef typename boost::graph_traits<G_t>::edge_descriptor Edge;

		typename boost::graph_traits<G_t>::vertex_iterator vIt, vEnd;
		for(boost::tie(vIt, vEnd) = boost::vertices(h); vIt != vEnd; ++vIt){
			typename boost::graph_traits<G_t>::adjacency_iterator nIt1, nIt2, nEnd;

//			typename boost::graph_traits<G_t>::in_edge_iterator in_nIt, in_nEnd;
//			for(boost::tie(in_nIt, in_nEnd) = boost::in_edges(v, _g); in_nIt != in_nEnd; ++in_nIt){

			for(boost::tie(nIt1, nEnd) = boost::adjacent_vertices(*vIt, h); nIt1 != nEnd; ++nIt1){
				for(boost::tie(nIt2, nEnd) = boost::adjacent_vertices(*vIt, h); nIt2 != nEnd; ++nIt2){
					if(nIt1 == nIt2){
						continue;
					}
					std::pair<Edge, bool> e = boost::edge(*nIt1, *nIt2, h);

					if(e.second){
						double w = boost::get(boost::edge_weight_t(), h, e.first);
						boost::add_edge(*nIt1, *nIt2, w, baseclass::g);
					}
					else{
						std::pair<Edge, bool> e1 = boost::edge(*nIt1, *vIt, h);
						if(e1.second){
							std::pair<Edge, bool> e2 = boost::edge(*vIt, *nIt2, h);
							double w1 = boost::get(boost::edge_weight_t(), h, e1.first);
							double w2 = boost::get(boost::edge_weight_t(), h, e2.first);

							boost::add_edge(*nIt1, *nIt2, abs(w1)+abs(w2), baseclass::g);
						}
						else{}
					}
				}
			}
		}

		typename boost::graph_traits<G_t>::edge_iterator eIt, eEnd;
		for(boost::tie(eIt, eEnd) = boost::edges(h); eIt != eEnd; ++eIt){
			std::pair<Edge, bool> e_g = boost::edge(boost::source(*eIt, h), boost::target(*eIt, h), baseclass::g);
			std::pair<Edge, bool> e_h = boost::edge(boost::source(*eIt, h), boost::target(*eIt, h), h);
			if(!e_g.second){
				double w = boost::get(boost::edge_weight_t(), h, e_h.first);
				boost::add_edge(boost::source(*eIt, h), boost::target(*eIt, h), w, baseclass::g);
			}
		}

		std::cout << "old edges: " << boost::num_edges(h) << std::endl;
		std::cout << "new edges: " << boost::num_edges(baseclass::g) << std::endl;

		baseclass::g = h;
	}
};

} //namespace

#endif //guard

/*
 * Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Author: Felix Salfelder, 2022
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

#ifndef UG4_GRAPH_INTERFACE_BIDIR_BOOST_H
#define UG4_GRAPH_INTERFACE_BIDIR_BOOST_H

#include "bidirectional.h"

namespace boost{

template <class T> struct graph_traits<ug::BidirectionalMatrix<ug::SparseMatrix<T>>>{
	typedef int vertex_descriptor;
	typedef SM_edge<T> edge_descriptor;
	typedef directed_tag bidirectional_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef SM_out_edge_iterator<T> out_edge_iterator;
	typedef SM_out_edge_iterator<T> in_edge_iterator;
//	typedef SM_out_edge_iterator<T> in_edge_iterator;
	typedef SM_adjacency_iterator<T> adjacency_iterator;
};

template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::BidirectionalMatrix<ug::SparseMatrix<T>> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<class T>
int out_degree(int v, ug::BidirectionalMatrix<T> const& M)
{
	return M.num_connections(v);
}

template<class T>
int in_degree(int v, ug::BidirectionalMatrix<T> const& M)
{
	return M.num_connections(v);
}

template<class T>
size_t source(SM_edge<T> const& e, ug::BidirectionalMatrix<ug::SparseMatrix<T>> const&)
{
	return e.row();
}

template<class T>
size_t target(SM_edge<T> const& e, ug::BidirectionalMatrix<ug::SparseMatrix<T>> const& m)
{
	return e.column(m);
}

template<class T>
std::pair<typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator,
          typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator>
adjacent_vertices(size_t v, ug::BidirectionalMatrix<ug::SparseMatrix<T>> const& M)
{
	typedef typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator a;

	typename ug::SparseMatrix<T>::const_row_iterator b = M.begin_row(v);
	typename ug::SparseMatrix<T>::const_row_iterator e = M.end_row(v);

	return std::make_pair(a(&b), a(&e));
}

template<class T>
std::pair<typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator,
          typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator>
coadjacent_vertices(size_t v, ug::BidirectionalMatrix<ug::SparseMatrix<T>> const& M)
{
	typedef typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator a;

	typename ug::SparseMatrix<T>::const_row_iterator b = M.begin_col(v);
	typename ug::SparseMatrix<T>::const_row_iterator e = M.end_col(v);

	return std::make_pair(a(&b), a(&e));
}


template<class T>
inline std::pair<SM_out_edge_iterator<T>, SM_out_edge_iterator<T>>
					out_edges(size_t v, ug::BidirectionalMatrix<ug::SparseMatrix<T>> const& g)
{
	typedef SM_out_edge_iterator<T> Iter;
   auto a = adjacent_vertices(v, g);
	return std::make_pair(Iter(v, a.first), Iter(v, a.second));
}

template<class T>
inline std::pair<SM_out_edge_iterator<T>, SM_out_edge_iterator<T>>
					in_edges(size_t v, ug::BidirectionalMatrix<ug::SparseMatrix<T>> const& g)
{
	typedef SM_out_edge_iterator<T> Iter;
   auto a = coadjacent_vertices(v, g);
	return std::make_pair(Iter(v, a.first), Iter(v, a.second));
}


} // boost

#endif // guard

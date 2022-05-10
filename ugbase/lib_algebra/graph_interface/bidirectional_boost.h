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

struct BS_traversal_tag
    : adjacency_graph_tag, bidirectional_graph_tag, vertex_list_graph_tag {};

template <class T>
struct graph_traits<ug::BidirectionalMatrix<T>>{
	typedef typename T::value_type value_type;
	typedef int vertex_descriptor;
	typedef SM_edge<value_type> edge_descriptor;
	typedef bidirectional_tag directed_category;
	typedef BS_traversal_tag traversal_category;
	typedef disallow_parallel_edge_tag edge_parallel_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef SM_out_edge_iterator<value_type> out_edge_iterator;
	typedef SM_out_edge_iterator<value_type> in_edge_iterator;
	typedef SM_adjacency_iterator<value_type> adjacency_iterator;
	typedef int degree_size_type;
	typedef int vertices_size_type;
};

#if 0
// deduplicate?
template <class T> struct graph_traits<ug::BidirectionalMatrix<ug::ParallelMatrix<ug::SparseMatrix<T>>>>{
	typedef int vertex_descriptor;
	typedef SM_edge<T> edge_descriptor;
	typedef bidirectional_tag directed_category;
	typedef BS_traversal_tag traversal_category;
	typedef disallow_parallel_edge_tag edge_parallel_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef SM_out_edge_iterator<T> out_edge_iterator;
	typedef SM_out_edge_iterator<T> in_edge_iterator;
	typedef SM_adjacency_iterator<T> adjacency_iterator;
	typedef int degree_size_type;
	typedef int vertices_size_type;
};
#endif

template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::BidirectionalMatrix<T> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<class T>
int num_vertices(ug::BidirectionalMatrix<T> const& M)
{
	return M.num_rows();
}

template<class T>
int out_degree(int v, ug::BidirectionalMatrix<T> const& M)
{
	return M.out_degree(v);
}

template<class T>
int in_degree(int v, ug::BidirectionalMatrix<T> const& M)
{
	return M.in_degree(v);
}

template<class T>
int degree(int v, ug::BidirectionalMatrix<T> const& M)
{ untested();
	return 2*M.degree(v);
}

template<class T>
size_t source(SM_edge<typename T::value_type> const& e, ug::BidirectionalMatrix<T> const&)
{
	return e.row();
}

template<class T>
size_t target(SM_edge<typename T::value_type> const& e, ug::BidirectionalMatrix<T> const& m)
{
	return e.column(m);
}

template<class T>
std::pair<typename graph_traits<T>::adjacency_iterator,
          typename graph_traits<T>::adjacency_iterator>
adjacent_vertices(size_t v, ug::BidirectionalMatrix<T> const& M)
{
	typedef typename graph_traits<T>::adjacency_iterator a;

	typename T::const_row_iterator b = M.begin_row(v);
	typename T::const_row_iterator e = M.end_row(v);

	return std::make_pair(a(&b, &e), a(&e, &e));
}

template<class T>
std::pair<typename graph_traits<T>::adjacency_iterator,
          typename graph_traits<T>::adjacency_iterator>
coadjacent_vertices(size_t v, ug::BidirectionalMatrix<T> const& M)
{
	typedef typename T::value_type value_type;
	typedef typename graph_traits<ug::SparseMatrix<value_type>>::adjacency_iterator a;

	typename T::const_row_iterator b = M.begin_col(v);
	typename T::const_row_iterator e = M.end_col(v);

	return std::make_pair(a(&b, &e), a(&e, &e));
}


template<class T>
inline std::pair<SM_out_edge_iterator<typename T::value_type>,
                 SM_out_edge_iterator<typename T::value_type>>
					out_edges(size_t v, ug::BidirectionalMatrix<T> const& g)
{
	typedef typename T::value_type value_type;
	typedef SM_out_edge_iterator<value_type> Iter;
   auto a = adjacent_vertices(v, g);
	return std::make_pair(Iter(a.first), Iter(a.second));
}

template<class T>
inline std::pair<SM_out_edge_iterator<typename T::value_type>,
                 SM_out_edge_iterator<typename T::value_type>>
					in_edges(size_t v, ug::BidirectionalMatrix<T> const& g)
{
	typedef typename T::value_type value_type;
	typedef SM_out_edge_iterator<value_type> Iter;
   auto a = coadjacent_vertices(v, g);
	return std::make_pair(Iter(a.first), Iter(a.second));
}

template<class T>
inline SM_edge_weight_map<typename T::value_type, ug::BidirectionalMatrix<T>>
get(edge_weight_t, ug::BidirectionalMatrix<T> const & g) {
	typedef typename T::value_type value_type;
	return SM_edge_weight_map<value_type, ug::BidirectionalMatrix<T>>(g);
}

template<class T>
struct property_map<ug::BidirectionalMatrix<ug::SparseMatrix<T>>, vertex_index_t>{
	typedef sparse_matrix_index_map<T> type;
	typedef type const_type;
};

template<class T>
inline typename property_map<ug::BidirectionalMatrix<ug::SparseMatrix<T>>, vertex_index_t>::const_type
get(vertex_index_t, ug::BidirectionalMatrix<T> const& m)
{
	return sparse_matrix_index_map<typename T::value_type>(m);
}

} // boost

#endif // guard

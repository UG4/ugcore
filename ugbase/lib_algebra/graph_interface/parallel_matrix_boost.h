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

#ifndef UG_PARALLEL_MATRIX_BOOST_H
#define UG_PARALLEL_MATRIX_BOOST_H
#ifdef UG_PARALLEL
#include "lib_algebra/parallelization/parallel_matrix.h"
#include "lib_algebra/graph_interface/parallel_matrix.h"
#include "lib_algebra/small_algebra/storage/fixed_array.h"


namespace boost{

template <class T> struct graph_traits<ug::BGLParallelMatrix<ug::ParallelMatrix<T>>>{
	typedef ug::BGLParallelMatrix<ug::ParallelMatrix<T>> type;
	typedef ug::detail::bglp_vertex_descriptor vertex_descriptor;
	typedef typename type::edge edge_descriptor;
	typedef directed_tag directed_category;
	typedef typename type::out_edge_iterator out_edge_iterator;
	typedef typename type::adjacency_iterator adjacency_iterator;
	typedef typename type::vertex_iterator vertex_iterator;
	typedef disallow_parallel_edge_tag edge_parallel_category;
	typedef SM_traversal_tag traversal_category;
	typedef int degree_size_type;
	typedef int vertices_size_type;
	//typedef typename ug::SparseMatrix<T>::const_row_iterator adjacency_iterator;
};

template<class T>
std::pair<typename ug::BGLParallelMatrix<T>::vertex_iterator,
          typename ug::BGLParallelMatrix<T>::vertex_iterator> vertices(
      ug::BGLParallelMatrix<T> const& M)
{
	auto b = M.begin_vertices();
	auto e = M.end_vertices();

	return std::make_pair(b,e);
}

#if 0
int owner(ug::detail::bglp_vertex_descriptor x)
{
	return x.owner();
}

int local(ug::detail::bglp_vertex_descriptor x)
{
	return x.local();
}
#else
// not needed.
using ug::detail::owner;
using ug::detail::local;
#endif

template<class T>
int out_degree(typename ug::BGLParallelMatrix<T>::vertex_descriptor v, ug::BGLParallelMatrix<T> const& M)
{
//	return M.out_degree(v);
	incomplete();
	return 17;
}

template<class T>
typename ug::BGLParallelMatrix<T>::vertex_descriptor
source(typename ug::BGLParallelMatrix<T>::edge const& e, ug::BGLParallelMatrix<T> const&)
{
	return e.source();
}

template<class T>
typename ug::BGLParallelMatrix<T>::vertex_descriptor
target(typename ug::BGLParallelMatrix<T>::edge const& e, ug::BGLParallelMatrix<T> const&)
{
	return e.target();
}

template<class T>
size_t num_vertices(ug::BGLParallelMatrix<T> const&)
{
	incomplete();
	return 0;
}

template<class T>
inline SM_edge_weight_map<T>
//inline typename property_map<ug::SparseMatrix<T>, edge_weight_t>::type
get(edge_weight_t, ug::BGLParallelMatrix<T> const & g) {
	incomplete();
//	return SM_edge_weight_map<T>(g);
}

template<class T>
inline std::pair<typename ug::BGLParallelMatrix<T>::out_edge_iterator,
                 typename ug::BGLParallelMatrix<T>::out_edge_iterator>
					out_edges(
							typename ug::BGLParallelMatrix<T>::vertex_descriptor v,
						  	ug::BGLParallelMatrix<T> const& M)
{
	auto l = local(v);
//	auto o = owner(v);

	auto b = M.begin_out_edges(l);
	auto e = M.end_out_edges(l);

	return std::make_pair(b, e);
}

template<class T>
inline std::pair<typename ug::BGLParallelMatrix<T>::adjacency_iterator,
                 typename ug::BGLParallelMatrix<T>::adjacency_iterator>
					adjacent_vertices(
							typename ug::BGLParallelMatrix<T>::vertex_descriptor v,
							ug::BGLParallelMatrix<T> const& M)
{
	auto l = local(v);
	auto b = M.begin_adjacent_vertices(l);
	auto e = M.end_adjacent_vertices(l);

	return std::make_pair(b,e);
}

template<class T>
class bglp_matrix_index_map
    : public put_get_helper<size_t, bglp_matrix_index_map<T> > { //
public:
	typedef size_t vertex_index_type;
	typedef size_t vertex_descriptor;
	typedef readable_property_map_tag category;
	typedef vertex_index_type value_type;
	typedef vertex_index_type reference;
	typedef vertex_descriptor key_type;

	bglp_matrix_index_map(bglp_matrix_index_map const& p) {
	}
	bglp_matrix_index_map(ug::SparseMatrix<T>const&, boost::vertex_index_t) { untested();
	}
	template<class X>
	bglp_matrix_index_map(X const&) {
	}
	template <class T_>
	value_type operator[](T_ x) const {
		return x.local();
	}
	value_type operator[](int x) const {
		return x;
	}
	bglp_matrix_index_map& operator=(const bglp_matrix_index_map& s) { untested();
		return *this;
	}
};

template<class T>
struct property_map<ug::BGLParallelMatrix<T>, vertex_index_t>{
	typedef bglp_matrix_index_map<T> type;
	typedef type const_type;
};

template<class T>
typename property_map<ug::BGLParallelMatrix<T>, vertex_index_t>::const_type
get(vertex_index_t, ug::BGLParallelMatrix<T> const& m)
{
	typedef typename property_map<ug::BGLParallelMatrix<T>, vertex_index_t>::type type;
	return type(m);
}

////////////////////////////////////////////////////////////////////////////////

// duplicate, same as graph_traits<ug::SparseMatrix<T>>?
template <class T> struct graph_traits<ug::ParallelMatrix<ug::SparseMatrix<T>>>{
	typedef int vertex_descriptor;
	typedef SM_edge<T> edge_descriptor;
	typedef directed_tag directed_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef SM_out_edge_iterator<T> out_edge_iterator;
	typedef SM_adjacency_iterator<T> adjacency_iterator;
	typedef disallow_parallel_edge_tag edge_parallel_category;
	typedef SM_traversal_tag traversal_category;
	typedef int degree_size_type;
	typedef int vertices_size_type;
	//typedef typename ug::SparseMatrix<T>::const_row_iterator adjacency_iterator;
};

template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::ParallelMatrix<ug::SparseMatrix<T>> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<class T>
struct property_map<ug::ParallelMatrix<T>, vertex_index_t>{
	typedef typename T::value_type value_type;
	typedef sparse_matrix_index_map<value_type> type;
	typedef type const_type;
};

} // boost

namespace ug {

using boost::counting_iterator;

#if 0
template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::ParallelMatrix<ug::SparseMatrix<T>> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<class T>
inline std::pair<boost::SM_out_edge_iterator<T>, boost::SM_out_edge_iterator<T>>
					out_edges(size_t v, ug::ParallelMatrix<ug::SparseMatrix<T>> const& g)
{
	typedef boost::SM_out_edge_iterator<T> Iter;
   auto a = boost::adjacent_vertices(v, g);
	return std::make_pair(Iter(a.first), Iter(a.second));
}
#endif

// stuff neeeded in boost::print_graph (?)
using boost::out_edges;
using boost::vertices;
using boost::source;
using boost::target;

} // ug

#endif // UG_PARALLEL
#endif // guard

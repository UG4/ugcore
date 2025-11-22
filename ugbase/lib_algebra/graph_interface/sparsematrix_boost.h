/* Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Felix Salfelder 2022
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 */

// BGL interface for cpu sparse matrix (dynamic CRS).

#ifndef UG_SPARSEMATRIX_BOOST_H
#define UG_SPARSEMATRIX_BOOST_H

#include "common/util/trace.h"
#include "lib_algebra/cpu_algebra/sparsematrix.h"

#include <boost/graph/properties.hpp> // put_get_helper
#include <boost/iterator/counting_iterator.hpp>

#include <cstdint> // for std::intmax_t

namespace boost{

using ug::iszero; // for now

// a bit of a pointer dance,
// because ug::SparseMatrix<T>::const_row_iterator is not default constructible.
template<typename T>
class SM_adjacency_iterator : public iterator_facade<
	SM_adjacency_iterator<T>,
	typename ug::SparseMatrix<T>::const_row_iterator,
	std::input_iterator_tag,
	size_t, // <= reference
	std::intmax_t // difference_type
	 >{ //
	using M = ug::SparseMatrix<T>;
	using iter_t = typename M::const_row_iterator;
	using base_class = iterator_facade<
		SM_adjacency_iterator,
		typename ug::SparseMatrix<T>::const_row_iterator,
		std::input_iterator_tag,
		size_t, // <= reference
		std::intmax_t // difference_type
	>;

public:
	using value_type = iter_t*;
	using difference_type = intmax_t;
	using reference = size_t;

public:
	SM_adjacency_iterator() : base_class(), _base(nullptr), _end(nullptr){
	}
	SM_adjacency_iterator(SM_adjacency_iterator&& p) = delete;
	SM_adjacency_iterator(SM_adjacency_iterator const& p)
		: _base(p._base?(new iter_t(*p._base)) : nullptr)
		, _end(p._end?(new iter_t(*p._end)) : nullptr){
	}
	SM_adjacency_iterator(value_type p, value_type e) : base_class(){
		assert(p);
		_base = new iter_t(*p);
		_end = new iter_t(*e);
		skip_zeroes();
	}
	~SM_adjacency_iterator() {
		delete _base;
		delete _end;
		_base = nullptr;
		_end = nullptr;
	}
	SM_adjacency_iterator& operator = (SM_adjacency_iterator&& other) = delete;
	SM_adjacency_iterator& operator = (const SM_adjacency_iterator& other) {
		if(other._base){
			delete _base;
			_base = new iter_t(*other._base);
		}else{ untested();
			_base = nullptr;
		}
		if(other._end){
			delete _end;
			_end = new iter_t(*other._end);
		}else{ untested();
			_end = nullptr;
		}
		return *this;
	}

	T const& value() const {
		assert(_base);
		return _base->value();
	}
	int row() const {
		assert(_base);
		return _base->row();
	}
	size_t idx() const { untested();
		assert(_base);
		return _base->idx();
	}
	reference dereference() const {
		assert(_base);
		return _base->index(); // row?
	}
	bool operator == (const SM_adjacency_iterator& other) const {
		assert(_base);
		assert(other._base);
		return *_base == *other._base;
	}
	bool operator != (const SM_adjacency_iterator& other) const {
		assert(_base);
		assert(other._base);
		return *_base != *other._base;
	}

private:
	// could use boost::filter_iterator, but does not work yet.
	void skip_zeroes(){
		while((*_base) != (*_end)){
			if(iszero(_base->value())){
				++(*_base);
			}else{
				break;
			}
		}
	}
	bool equal(SM_adjacency_iterator const& other) const { untested();
		assert(_base);
		assert(other._base);
		return *_base == *other._base;
	}
	void increment() {
		assert(_base);
		assert(_end);
		++(*_base);
		skip_zeroes();
	}
	void decrement() { untested();
		// incomplete(); // don't use. too complicated.
		assert(_base);
		--(*_base);
	}

private:
	value_type _base;
	value_type _end;
	friend class iterator_core_access;
}; // SM_adjacency_iterator

template<typename T>
class SM_edge{
public:
//	SM_edge(size_t v, size_t w) : _row(v), _idx(w) { untested();
//	}
	SM_edge(SM_edge const&) = default;
	SM_edge(SM_edge&&) = default;
	explicit SM_edge() { untested();
	}
	SM_edge(int s, int t, T const& v) : _source(s), _target(t), _value(v) {
	}
	SM_edge& operator = (SM_edge const&) = default;
	SM_edge& operator = (SM_edge&&) = default;

public:
	int source() const{
		return _source;
	}
	int target() const{
		return _target;
	}
	template <typename X>
	T const& value(X const&) const{
		return _value;
	}

public:
	bool operator == (const SM_edge& other) const { untested();
		return _source == other._source && _target == other._target;
	}
	bool operator != (const SM_edge& other) const { untested();
		return !operator == (other);
	}

private:
	int _source;
	int _target;
	T _value;
}; // SM_edge

template<typename T, typename M>
int source(SM_edge<T> const& e, M const&)
{ untested();
	return e.source();
}

template<typename T, typename M>
int target(SM_edge<T> const& e, M const&)
{ untested();
	return e.target();
}

template<typename T, bool out=true>
class SM_out_edge_iterator : public iterator_facade<
	SM_out_edge_iterator<T, out>,
	SM_adjacency_iterator<T>,
	// bidirectional_traversal_tag, // breaks InputIterator (why?)
	std::input_iterator_tag,
	SM_edge<T>, // <= reference
	std::intmax_t // difference_type
	 >{ //
public: // types
	using base_class = iterator_facade<
		SM_out_edge_iterator<T, out>,
		SM_adjacency_iterator<T>,
		std::input_iterator_tag,
		SM_edge<T>, // <= reference
		std::intmax_t // difference_type
	>;
	using M = ug::SparseMatrix<T>;
	using value_type = SM_adjacency_iterator<T>;
	using difference_type = intmax_t;
	using reference = SM_edge<T>;
	using edge_type = SM_edge<T>;

public: // construct
	explicit SM_out_edge_iterator() : base_class() {
	}
	explicit SM_out_edge_iterator(SM_adjacency_iterator<T> w, int src)
	    : base_class(), _base(w), _src(src) {
		if(out){
		}else{
		}
	}
	/* explicit */ SM_out_edge_iterator(SM_out_edge_iterator const& p)
	    : base_class(p), _base(p._base), _src(p._src){
		if(out){
		}else{
		}
	}
	SM_out_edge_iterator(SM_out_edge_iterator&& p) = delete; // why?
	~SM_out_edge_iterator() = default;

	SM_out_edge_iterator& operator = (SM_out_edge_iterator&& other) =  default;
	SM_out_edge_iterator& operator = (SM_out_edge_iterator const& other) = default;
#if 0
public: // op
	SM_out_edge_iterator operator + (int a) const{ untested();
		SM_out_edge_iterator ret(*this);
		ret.base.first += a;
		return ret;
	}
	SM_out_edge_iterator operator - (int a) const{ untested();
		SM_out_edge_iterator ret(*this);
		ret.base.first -= a;
		return ret;
	}
	difference_type operator - (SM_out_edge_iterator const& other) const{ untested();
		SM_out_edge_iterator ret(*this);
		return ret.base.first - other.base.first;
	}
#endif
private:
	reference dereference() const {
		if(out){
			return edge_type(_src, *_base, _base.value());
		}else{
			return edge_type(*_base, _src, _base.value());
		}
	}
	bool equal(SM_out_edge_iterator const& other) const { untested();
		assert(_base.first == other._base.first);
		return _base.second == other._base.second;
	}
	void increment() { untested();
		++_base;
	}
	void decrement() { untested();
		--_base;
	}
	//			bool operator == (const SM_out_edge_iterator& other) const
	//			{ incomplete();
	//				return false;
	//			}
public:
	SM_out_edge_iterator& operator ++ () {
		++_base;
		return *this;
	}
	SM_out_edge_iterator operator ++ (int) { untested();
		SM_out_edge_iterator copy(*this);
		++*this;
		return copy;
	}
	bool operator == (const SM_out_edge_iterator& other) const {
		return _base == other._base;
	}
	bool operator != (const SM_out_edge_iterator& other) const {
		return _base != other._base;
	}

private:
	value_type _base;
	int _src;
	friend class iterator_core_access;
}; // SM_out_edge_iterator

struct SM_traversal_tag
    : adjacency_graph_tag, bidirectional_graph_tag, vertex_list_graph_tag {};

template <typename T> struct graph_traits<ug::SparseMatrix<T>>{
	using vertex_descriptor = int;
	using edge_descriptor = SM_edge<T>;
	using directed_category = directed_tag;
	using edge_parallel_category = disallow_parallel_edge_tag;
	using traversal_category = SM_traversal_tag;
	using vertex_iterator = counting_iterator<size_t>;
	using out_edge_iterator = SM_out_edge_iterator<T>;
	using adjacency_iterator = SM_adjacency_iterator<T>;
	using degree_size_type = int;
	using vertices_size_type = int;
};

template <typename T>
std::pair<typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator,
          typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator>
inline adjacent_vertices(size_t v, ug::SparseMatrix<T> const& M)
{
	assert(v<M.num_rows());
	using a = typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator;

	typename ug::SparseMatrix<T>::const_row_iterator b(M.begin_row(v));
	typename ug::SparseMatrix<T>::const_row_iterator e(M.end_row(v));

	return std::make_pair(a(&b, &e), a(&e, &e));
}

template <typename T>
inline std::pair<SM_out_edge_iterator<T>, SM_out_edge_iterator<T>>
					out_edges(int v, ug::SparseMatrix<T> const& g)
{
	assert(size_t(v)<g.num_rows());
	using Iter = SM_out_edge_iterator<T>;
	auto a = adjacent_vertices(v, g);
	return std::make_pair(Iter(a.first, v), Iter(a.second, v));
}

template<typename T>
class sparse_matrix_index_map
    : public put_get_helper<size_t, sparse_matrix_index_map<T> > { //
public:
	using vertex_index_type = size_t;
	using vertex_descriptor = size_t;
	using category = readable_property_map_tag;
	using value_type = vertex_index_type;
	using reference = vertex_index_type;
	using key_type = vertex_descriptor;

	sparse_matrix_index_map(sparse_matrix_index_map const& p) {
	}
	sparse_matrix_index_map(ug::SparseMatrix<T>const&, vertex_index_t) { untested();
	}
	template<typename X>
	sparse_matrix_index_map(X const&) {
	}
	template <typename T_>
	value_type operator [] (T_ x) const {
		return x;
	}
	sparse_matrix_index_map& operator = (const sparse_matrix_index_map& s) { untested();
		return *this;
	}
};

template<typename T>
struct property_map<ug::SparseMatrix<T>, vertex_index_t>{
	using type = sparse_matrix_index_map<T>;
	using const_type = type;
};

template <typename T>
typename property_map<ug::SparseMatrix<T>, vertex_index_t>::const_type
get(vertex_index_t, ug::SparseMatrix<T> const& m)
{
	return sparse_matrix_index_map<T>(m);
}

template<typename T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::SparseMatrix<T> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(M.num_rows());

	return std::make_pair(b,e);
}

template<typename T>
int num_vertices(ug::SparseMatrix<T> const& M)
{
	assert(M.num_rows() == M.num_cols());
	return M.num_rows();
}

template<typename T>
int out_degree(int v, ug::SparseMatrix<T> const& M)
{
	int c = 0;
	auto i = out_edges(v, M);
	for(; i.first != i.second; ++i.first) {
		++c;
	}
	return c;
}

// template<typename T>
// T get(edge_weight_t, ug::SparseMatrix<T> const& M, SM_edge<T> e)
// { untested();
// 	return e.value();
// }

template<typename T, typename M=ug::SparseMatrix<T>>
class SM_edge_weight_map :
	public put_get_helper< T /*algebraic connection?*/, SM_edge_weight_map<T> > { //
public:
	using edge_weight_type = T;
	using vertex_descriptor = int;
	using category = readable_property_map_tag;
	using value_type = edge_weight_type;
	using reference = T&;
	using key_type = vertex_descriptor;

	SM_edge_weight_map(SM_edge_weight_map const& p) : _g(p._g) { untested();
	}
	SM_edge_weight_map(ug::SparseMatrix<T>const & g, boost::edge_weight_t) : _g(g) { untested();
	}
	// bug?
	SM_edge_weight_map(M const& g) : _g(g) {
	}
	template <typename X>
	value_type operator [] (X const& x) const {
		//				assert(x == _g.position(x));
		return x.value(_g);
	}
	SM_edge_weight_map& operator = (const SM_edge_weight_map& s) { untested();
		assert(&s._g==&_g); (void)s;
		return *this;
	}
private:
	M const& _g;
}; // SM_edge_weight_map

// g++ 'enumeral_type' in template unification not implemented workaround
template<typename T, typename Tag>
struct property_map<ug::SparseMatrix<T>, Tag> {
//	using map_gen = typename gala_graph_property_map<Tag>::template bind_<T>;
//	using type = typename map_gen::type;
//	using const_type = typename map_gen::const_type;
};

template<typename T>
inline SM_edge_weight_map<T>
//inline typename property_map<ug::SparseMatrix<T>, edge_weight_t>::type
get(edge_weight_t, ug::SparseMatrix<T> const & g) {
	return SM_edge_weight_map<T>(g);
}

#if 0
template<typename T>
inline typename property_map<ug::SparseMatrix<T>, edge_all_t>::type
get(edge_all_t, ug::SparseMatrix<T> & g) { incomplete();
	using pmap_type = typename property_map<ug::SparseMatrix<T>, edge_all_t>::type;
	return pmap_type(&g);
}
#endif

} // boost

namespace ug{

// used from boost::print_graph, graph_utility.hpp.
// must be in ug, because of ADL. why don't they call boost::{out_edges,vertices}?
using boost::vertices;
using boost::out_edges;

}// ug

#endif

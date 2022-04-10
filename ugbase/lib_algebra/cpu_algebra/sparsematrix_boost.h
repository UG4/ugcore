/* Copyright (c) 2022:  G-CSC, Goethe University Frankfurt
 * Felix Salfelder 2022
 * This file is part of UG4.
 *
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 */

// BGL interface for cpu sparse martrix (dynamic CRS).

#include "sparsematrix.h"
#include <boost/graph/properties.hpp> // put_get_helper
#include <boost/iterator/counting_iterator.hpp>

namespace boost{

// a bit of a pointer dance,
// because ug::SparseMatrix<T>::const_row_iterator is not default constructible.
template<class T>
class SM_adjacency_iterator : public iterator_facade<
	SM_adjacency_iterator<T>,
	typename ug::SparseMatrix<T>::const_row_iterator,
	std::input_iterator_tag,
	size_t, // <= reference
	std::intmax_t // difference_type
	 >{
	typedef ug::SparseMatrix<T> M;
	typedef typename M::const_row_iterator iter_t;

public:
	typedef iter_t* value_type;
	typedef intmax_t difference_type;
	typedef size_t reference;

public:
	SM_adjacency_iterator() : _base(nullptr){}
	SM_adjacency_iterator(SM_adjacency_iterator&& p) = delete;
	SM_adjacency_iterator(SM_adjacency_iterator const& p)
		: _base(p._base?(new iter_t(*p._base)) : nullptr){
	}
	SM_adjacency_iterator(value_type p){
		assert(p);
		_base = new iter_t(*p);
	}
	~SM_adjacency_iterator(){
		delete _base;
		_base = nullptr;
	}
	size_t idx() const {
		assert(_base);
		return _base->idx();
	}
	reference dereference() const {
		assert(_base);
		return _base->index(); // row?
	}
	bool operator!=(const SM_adjacency_iterator& other) const {
		assert(_base);
		assert(other._base);
		return *_base != *other._base;
	}
	SM_adjacency_iterator operator=(const SM_adjacency_iterator& other) {
		if(other._base){
			_base = new iter_t(*other._base);
		}else{
			_base = nullptr;
		}
		return *this;
	}

private:
	bool equal(SM_adjacency_iterator const& other) const {
		assert(_base);
		assert(other._base);
		return *_base == *other._base;
	}
	void increment() {
		assert(_base);
		++(*_base);
	}
	void decrement() {
		assert(_base);
		--(*_base);
	}

private:
	value_type _base;
	friend class iterator_core_access;
}; // SM_adjacency_iterator

template<class T>
class SM_edge{
public:
	SM_edge(size_t v, size_t w) : _row(v), _idx(w) {}

public:
	size_t row() const{
		return _row;
	}
	size_t column(ug::SparseMatrix<T> const& m) const{
		return m.col(_idx);
	}

private:
	int _row;
	size_t _idx;
};

template<class T>
size_t source(SM_edge<T> const& e, ug::SparseMatrix<T> const&)
{
	return e.row();
}

template<class T>
size_t target(SM_edge<T> const& e, ug::SparseMatrix<T> const& m)
{
	return e.column(m);
}

template<class T>
class SM_out_edge_iterator : public iterator_facade<
	SM_out_edge_iterator<T>,
	std::pair<size_t, SM_adjacency_iterator<T>>,
	// bidirectional_traversal_tag, // breaks InputIterator (why?)
	std::input_iterator_tag,
	SM_edge<T>, // <= reference
	std::intmax_t // difference_type
	 >{
public: // types
	typedef ug::SparseMatrix<T> M;
	typedef typename std::pair<size_t, SM_adjacency_iterator<T> > value_type;
	typedef intmax_t difference_type;
	typedef SM_edge<T> reference;
	typedef SM_edge<T> edge_type;

public: // construct
	explicit SM_out_edge_iterator() {
	}
	explicit SM_out_edge_iterator(size_t v, SM_adjacency_iterator<T> w)
	    : _base(v, w) {
	}
	explicit SM_out_edge_iterator(SM_out_edge_iterator const& p)
	    : _base(p._base){
	}
	~SM_out_edge_iterator(){
	}
#if 0
public: // op
	SM_out_edge_iterator operator+(int a) const{ untested();
		SM_out_edge_iterator ret(*this);
		ret.base.first += a;
		return ret;
	}
	SM_out_edge_iterator operator-(int a) const{ untested();
		SM_out_edge_iterator ret(*this);
		ret.base.first -= a;
		return ret;
	}
	difference_type operator-(SM_out_edge_iterator const& other) const{ untested();
		SM_out_edge_iterator ret(*this);
		return ret.base.first - other.base.first;
	}
#endif
private:
	reference dereference() const {
		return edge_type(_base.first, _base.second.idx());
	}
	bool equal(SM_out_edge_iterator const& other) const {
		assert(_base.first == other._base.first);
		return _base.second == other._base.second;
	}
	void increment() {
		++_base.second;
	}
	void decrement() {
		--_base.second;
	}
	//			bool operator==(const SM_out_edge_iterator& other) const
	//			{ incomplete();
	//				return false;
	//			}
public:
	bool operator==(const SM_out_edge_iterator& other) const {
		return _base.second == other._base.second;
	}
	bool operator!=(const SM_out_edge_iterator& other) const {
		return _base.second != other._base.second;
	}

private:
	value_type _base;
	friend class iterator_core_access;
}; // SM_out_edge_iterator

template <class T> struct graph_traits<ug::SparseMatrix<T>>{
	typedef size_t vertex_descriptor;
	typedef SM_edge<T> edge_descriptor;
	typedef directed_tag directed_category;
	typedef counting_iterator<size_t> vertex_iterator;
	typedef SM_out_edge_iterator<T> out_edge_iterator;
	typedef SM_adjacency_iterator<T> adjacency_iterator;
	//typedef typename ug::SparseMatrix<T>::const_row_iterator adjacency_iterator;
};

template<class T>
std::pair<typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator,
          typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator>
adjacent_vertices(size_t v, ug::SparseMatrix<T> const& M)
{
	typedef typename graph_traits<ug::SparseMatrix<T>>::adjacency_iterator a;

	typename ug::SparseMatrix<T>::const_row_iterator b = M.begin_row(v);
	typename ug::SparseMatrix<T>::const_row_iterator e = M.end_row(v);

	return std::make_pair(a(&b), a(&e));
}

template<class T>
inline std::pair<SM_out_edge_iterator<T>, SM_out_edge_iterator<T>>
					out_edges(size_t v, ug::SparseMatrix<T> const& g)
{
	typedef SM_out_edge_iterator<T> Iter;
//	ug::SparseMatrix<T>* G = const_cast<ug::SparseMatrix<T>*>(&g); // HACK
   auto a = adjacent_vertices(v, g);
	return std::make_pair(Iter(v, a.first), Iter(v, a.second));
}

template<class T>
class sparse_matrix_index_map
	 : public put_get_helper<size_t, sparse_matrix_index_map<T> > {
	public:
		typedef size_t vertex_index_type;
		typedef size_t vertex_descriptor;
		typedef readable_property_map_tag category;
		typedef vertex_index_type value_type;
		typedef vertex_index_type reference;
		typedef vertex_descriptor key_type;

		sparse_matrix_index_map(sparse_matrix_index_map const& p) { }
		sparse_matrix_index_map(ug::SparseMatrix<T>const&, boost::vertex_index_t) { }
		sparse_matrix_index_map(ug::SparseMatrix<T>const&) { }
		template <class T_>
		value_type operator[](T_ x) const {
			return x;
		}
		sparse_matrix_index_map& operator=(const sparse_matrix_index_map& s) {
			return *this;
		}
};

template<class T>
struct property_map<ug::SparseMatrix<T>, vertex_index_t>{
	typedef sparse_matrix_index_map<T> type;
	typedef type const_type;
};

template<class T>
inline typename property_map<ug::SparseMatrix<T>, vertex_index_t>::const_type
get(vertex_index_t, ug::SparseMatrix<T> const& m){
	return sparse_matrix_index_map<T>(m);
}

template<class T>
std::pair<counting_iterator<size_t>, counting_iterator<size_t> > vertices(
      ug::SparseMatrix<T> const& M)
{
	counting_iterator<size_t> b(0);
	counting_iterator<size_t> e(7);

	return std::make_pair(b,e);
}

} // boost

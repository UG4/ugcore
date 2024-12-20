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

#ifndef UG_GRAPH_INTERFACE_PARALLEL_MATRIX_H
#define UG_GRAPH_INTERFACE_PARALLEL_MATRIX_H

#include "sparsematrix_boost.h"
#include <boost/iterator/filter_iterator.hpp>

namespace ug{

namespace detail{

class bglp_vertex_descriptor : public std::pair<int, int>{
public:
	bglp_vertex_descriptor() {}
	bglp_vertex_descriptor(int a, int b) : std::pair<int, int>(a, b) {
	}
public:
	int owner() const{ return first; }
	int local() const{ return second; }
};

inline int owner(bglp_vertex_descriptor p)
{
	return p.owner();
}
inline int local(bglp_vertex_descriptor p)
{
	return p.local();
}

} // detail

// BGL style access to matrix
// assign vertices to processes,
// provide access to cross process entries
template<class T>
class BGLParallelMatrix {
public: // types
	typedef detail::bglp_vertex_descriptor vertex_descriptor;
	typedef typename T::const_row_iterator const_row_iterator;
private:
	typedef std::vector<int> owners;
	typedef std::vector<int> ghosts;
	typedef typename boost::graph_traits<T>::vertex_iterator base_vertex_iterator;
	typedef typename boost::graph_traits<T>::out_edge_iterator base_edge_iterator;
	class vertex_iterator_ // facade?
	{ //

	public:

		using iterator_category = std::input_iterator_tag;
		using value_type = vertex_descriptor;
		using difference_type = ptrdiff_t;
		using pointer = vertex_descriptor*;
		using reference = vertex_descriptor;

		explicit vertex_iterator_() : _owners(nullptr) {}
		explicit vertex_iterator_(base_vertex_iterator b, owners const* o) : _base(b), _owners(o) {
		}
		bool operator==(vertex_iterator_ const& o) const{
			return _base == o._base;
		}
		bool operator!=(vertex_iterator_ const& p) const{
			return !operator==(p);
		}
		vertex_iterator_& operator++(int) {
			incomplete();
			return *this;
		}
		vertex_iterator_& operator++() {
			++_base;
			return *this;
		}
		vertex_descriptor operator*() const {
			auto b = *_base;
			assert(_owners);
			return vertex_descriptor((*_owners)[b], *_base);
		}
	private:
		base_vertex_iterator _base;
		owners const* _owners;
		ghosts const* _ghosts;
	};

	struct filter_local{
		bool operator()(vertex_descriptor v) const{
			int rank;
#ifdef UG_PARALLEL
			MPI_Comm_rank(MPI_COMM_WORLD, &rank);
			return v.owner() == rank;
#else
			return true;
#endif
		}
	};

public:
	typedef boost::filter_iterator<filter_local, vertex_iterator_> vertex_iterator;
	class adjacency_iterator // facade?
		/* todo remove deprecated inheritance
		 :public std::iterator<std::input_iterator_tag, vertex_descriptor, ptrdiff_t, vertex_descriptor, vertex_descriptor> */
	{ //
	public:

		using iterator_category = std::input_iterator_tag;
		using value_type = vertex_descriptor;
		using difference_type = ptrdiff_t;
		using pointer = vertex_descriptor*;
		using reference = vertex_descriptor;

		typedef typename boost::graph_traits<T>::adjacency_iterator base;
	public:
		adjacency_iterator()
		    : _owners(nullptr), _ghosts(nullptr) {
		}
		adjacency_iterator(base b, owners const* o, ghosts const* g)
		    : _base(b), _owners(o), _ghosts(g) {
		}
	public:
		bool operator==(adjacency_iterator const& o) const{
			return _base == o._base;;
		}
		bool operator!=(adjacency_iterator const& p) const{
			return !operator==(p);
		}
		adjacency_iterator& operator++(int) {
			incomplete();
			return *this;
		}
		adjacency_iterator& operator++() {
			++_base;
			return *this;
		}
		vertex_descriptor operator*() const {
			int i = *_base;
			int o = (*_owners)[i];
			int l = (*_ghosts)[i];

			return vertex_descriptor(o, l);
		}
		base _base;
		owners const* _owners;
		ghosts const* _ghosts;
	}; // adjacency_iterator

	class edge;
	class out_edge_iterator : public boost::iterator_facade<
	    out_edge_iterator,
	    base_edge_iterator,
		 std::input_iterator_tag,
		 edge, // <= reference
		 std::intmax_t // difference_type
	 >{ //
	public:
		out_edge_iterator()
		    : _owners(nullptr), _ghosts(nullptr), _matrix(nullptr) {
		}
		out_edge_iterator(base_edge_iterator b, owners const* o, ghosts const* g, T const* m)
		    : _base(b), _owners(o), _ghosts(g), _matrix(m) {
		}
		bool operator==(out_edge_iterator const& p) const{
			return _base == p._base;
		}
		bool operator!=(out_edge_iterator const& p) const{
			return !operator==(p);
		}
		out_edge_iterator& operator++() {
			++_base;
			return *this;
		}
		edge operator*() const;
	private:
		base_edge_iterator _base;
		owners const* _owners;
		ghosts const* _ghosts;
		T const* _matrix;
	};
	class edge {
	public:
		edge(vertex_descriptor s, vertex_descriptor t) : _s(s), _t(t) {
		}
	public:

		vertex_descriptor source() const{
			return _s;
		}
		vertex_descriptor target() const{
			return _t;
		}

	private:
		vertex_descriptor _s;
		vertex_descriptor _t;

	}; // edge

public:
	explicit BGLParallelMatrix(T const* m=nullptr)
	    : _matrix(m) {
		if(m){
			refresh();
		}else{ untested();
		}
	}
//	explicit BGLParallelMatrix(ug::ParallelMatrix<T> const& o)
//	    : _matrix(o._matrix), _matrix_transpose(o._matrix_transpose) { untested();
//	}
	BGLParallelMatrix& operator=(BGLParallelMatrix const& o) { untested();
		_matrix = o._matrix;
		_matrix_transpose = o._matrix_transpose;
		return *this;
	}

public: // interface
	void refresh();

	int num_rows() const {
		return std::max(_matrix->num_rows(), _matrix->num_cols());
	}
	int num_cols() const { untested();
		return num_rows();
	}
	int num_connections(int v) const {
		if(v<_matrix->num_rows()){
			return _matrix->num_connections(v);
		}else{
			assert(v<_matrix_transpose.num_rows());
			return _matrix_transpose.num_connections(v);
		}
	}
	int out_degree(int v) const {
		// could call in_degree?
		return in_degree(v);
	}
	int in_degree(int v) const {
		// could use difference, requires zero-pruning in _matrix_transpose.
		if(v<_matrix_transpose.num_rows()){
			return boost::out_degree(v, _matrix_transpose);
		}else{
			assert(_matrix);
			assert(v<_matrix->num_rows());
			return boost::out_degree(v, *_matrix);
		}
	}
	int degree(int v) const { untested();
		return 2*out_degree(v);
	}

	adjacency_iterator begin_adjacent_vertices(int row) const { untested();
		assert(_matrix);
		auto p = boost::adjacent_vertices(row, *_matrix);
		return adjacency_iterator(p.first, &_owners, &_ghosts);
	}
	adjacency_iterator end_adjacent_vertices(int row) const { untested();
		assert(_matrix);
		auto p = boost::adjacent_vertices(row, *_matrix);
		return adjacency_iterator(p.second, &_owners, &_ghosts);
	}
	out_edge_iterator begin_out_edges(int row) const {
		assert(_matrix);
		auto p = boost::out_edges(row, *_matrix);
		return out_edge_iterator(p.first, &_owners, &_ghosts, _matrix);
	}
	out_edge_iterator end_out_edges(int row) const {
		assert(_matrix);
		auto p = boost::out_edges(row, *_matrix);
		return out_edge_iterator(p.second, &_owners, &_ghosts, _matrix);
	}

#if 0
	const_row_iterator begin_row(int row) const {
		assert(_matrix);
		if(row<_matrix->num_rows()){
		}else{
			row = 0;
		}
		return _matrix->begin_row(row);
	}
	const_row_iterator end_row(int row) const {
		assert(_matrix);
		if(row<_matrix->num_rows()){
			return _matrix->end_row(row);
		}else{
			return _matrix->begin_row(0);
		}
	}
#endif

	vertex_iterator begin_vertices() const {
		assert(_matrix);
		auto v = boost::vertices(*_matrix);
		auto b = vertex_iterator_(v.first, &_owners);
		auto e = vertex_iterator_(v.second, &_owners);
		// filter_local f(_owners);
		return vertex_iterator( b, e );
	}

	vertex_iterator end_vertices() const {
		assert(_matrix);
		auto v = boost::vertices(*_matrix);
		auto e = vertex_iterator_(v.second, &_owners);
		return vertex_iterator( e, e );
	}

	const_row_iterator begin_col(int col) const {
		if(col<_matrix_transpose.num_rows()){
		}else{
			col = 0;
		}
		return _matrix_transpose.begin_row(col);
	}
	const_row_iterator end_col(int col) const {
		if(col<_matrix_transpose.num_rows()){
			return _matrix_transpose.end_row(col);
		}else{
			return _matrix_transpose.begin_row(0);
		}
	}

private:
public: // BUG
	T const* _matrix;
	owners _owners;
	ghosts _ghosts;
private:
	T _matrix_transpose;
}; // BGLParallelMatrix

template<class T>
typename BGLParallelMatrix<T>::edge BGLParallelMatrix<T>::out_edge_iterator::operator*() const
{
	auto e = *_base;
	assert(_matrix);
	auto s = boost::source(e, _matrix);
	auto t = boost::target(e, _matrix);
	typedef typename BGLParallelMatrix<T>::edge E;
	vertex_descriptor v((*_owners)[s], (*_ghosts)[s]);
	vertex_descriptor w((*_owners)[t], (*_ghosts)[t]);

	return E(v, w);
}

template<class T>
void BGLParallelMatrix<T>::refresh()
{
	assert(_matrix);
#ifdef UG_PARALLEL
	assert(_matrix->num_rows() == _matrix->num_cols()); // for now.
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	_owners.resize(_matrix->num_rows(), rank);
	_ghosts.resize(_matrix->num_rows()); // map slave nodes to master local index

	int l = 0;
	for(auto& i : _ghosts){
		i = l++;
	}

	auto L = _matrix->layouts();
	assert(L);

	pcl::ProcessCommunicator pc = L->proc_comm();
	// void PC::send_data(void* pBuffer, int bufferSize, int destProc, int tag) const;

	// send indexmap to slaves
	IndexLayout const& idx_master = L->master();
	for(auto k = idx_master.begin(); k!=idx_master.end(); ++k){
		IndexLayout::Interface const& I = idx_master.interface(k);
		std::vector<int> mmap(I.size());
		int dest = idx_master.proc_id(k);
		auto b = mmap.begin();

		for(auto i : I){
			*b = i.elem;
			++b;
		}
		pc.send_data(mmap.data(), mmap.size()*sizeof(int), dest, 17);
	}

	IndexLayout const& idx_slave = L->slave();
	for(auto k = idx_slave.begin(); k!=idx_slave.end(); ++k){
		IndexLayout::Interface const& I = idx_slave.interface(k);
		int srank = idx_slave.proc_id(k);

		std::vector<int> mmap(I.size());
		pc.receive_data(mmap.data(), mmap.size()*sizeof(int), srank, 17);
//		for(auto i : mmap){
//			std::cerr << "received " << rank << " from " << srank << " " << i << "\n";
//		}

		int j=0;
		for(auto i : I){
			assert(_owners[i.elem] == rank); // only one owner per node.
			_owners[i.elem] = srank;
			++j;
			_ghosts[i.elem] = mmap[j];
		}
	}
#endif
}

} // ug

#endif // guard

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

#ifndef UG_GRAPH_INTERFACE_BIDIR_H
#define UG_GRAPH_INTERFACE_BIDIR_H

#include "sparsematrix_boost.h"

#ifndef NDEBUG
#include "common/stopwatch.h"
#endif

namespace ug{

// give access to edges both ways.
// keep a transpose matrix to cache iterators.
template<class T>
class BidirectionalMatrix{
public: // types
	typedef typename T::const_row_iterator const_row_iterator;
public:
	explicit BidirectionalMatrix(T const* m=nullptr)
	    : _matrix(m) {
		if(m){
			refresh();
		}else{
		}
	}
	explicit BidirectionalMatrix(BidirectionalMatrix const& o)
	    : _matrix(o._matrix) {
		_matrix_transpose = o._matrix_transpose; // BUG: missing copy constructor.
	}
	BidirectionalMatrix& operator=(BidirectionalMatrix const& o) {
		_matrix = o._matrix;
		_matrix_transpose = o._matrix_transpose;
		return *this;
	}

public: // interface
	void refresh(){
#ifndef NDEBUG
		double t = -get_clock_s();
#endif
		assert(_matrix);
		_matrix_transpose.set_as_transpose_of2(*_matrix);
		assert(_matrix->num_rows() == _matrix_transpose.num_cols());
		assert(_matrix->num_cols() == _matrix_transpose.num_rows());
#ifndef NDEBUG
		t += get_clock_s();
		UG_LOG("Transpose in bidirectional refresh " << t << "\n");
#endif
	}

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
		assert(_matrix);
		if(size_t(v)<_matrix->num_rows()){
			assert(size_t(v)<_matrix->num_rows());
			return boost::out_degree(v, *_matrix);
		}else{
			return 0;
		}
	}
	int in_degree(int v) const {
		// could use difference, requires zero-pruning in _matrix_transpose.
		if(size_t(v)<_matrix_transpose.num_rows()){
			return boost::out_degree(v, _matrix_transpose);
		}else{
			return 0;
		}
	}
	int degree(int v) const { untested();
		return 2*out_degree(v);
	}

	const_row_iterator begin_row(int row) const {
		assert(_matrix);
		if(size_t(row)<_matrix->num_rows()){
		}else{
			row = 0;
		}
		return _matrix->begin_row(row);
	}
	const_row_iterator end_row(int row) const {
		assert(_matrix);
		if(size_t(row)<_matrix->num_rows()){
			return _matrix->end_row(row);
		}else{
			return _matrix->begin_row(0);
		}
	}

	const_row_iterator begin_col(int col) const {
		if(size_t(col)<_matrix_transpose.num_rows()){
		}else{
			col = 0;
		}
		return _matrix_transpose.begin_row(col);
	}
	const_row_iterator end_col(int col) const {
		if(size_t(col)<_matrix_transpose.num_rows()){
			return _matrix_transpose.end_row(col);
		}else{
			return _matrix_transpose.begin_row(0);
		}
	}

private:
	T const* _matrix;
	T _matrix_transpose;
};

} // ug

#endif // guard

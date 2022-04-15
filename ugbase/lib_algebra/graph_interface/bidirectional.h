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

#ifndef UG4_GRAPH_INTERFACE_BIDIR_H
#define UG4_GRAPH_INTERFACE_BIDIR_H

#include "trace.h"

namespace ug{

// give access to edges both ways.
// keep a transpose matrix to cache iterators.
template<class T>
class BidirectionalMatrix{
public: // types
	typedef typename T::const_row_iterator const_row_iterator;
public:
	explicit BidirectionalMatrix(ConstSmartPtr<T> m) : _matrix(m) { untested();
		refresh();
	}
	explicit BidirectionalMatrix(T const* m) : _matrix(nullptr) {
		T* m_ = new T;
		m_->set_as_copy_of(*m);
		_matrix = make_sp(m_);
		refresh();
	}

public: // interface
	void refresh(){
		assert(_matrix);
		_matrix_transpose.set_as_transpose_of(*_matrix);
		assert(_matrix->num_rows() == _matrix_transpose.num_cols());
		assert(_matrix->num_cols() == _matrix_transpose.num_rows());
	}

	int num_rows() const {
		return _matrix->num_rows();
	}
	int num_cols() const { untested();
		return _matrix->num_cols();
	}
	int num_connections(int v) const {
		return _matrix->num_connections(v);
	}

	const_row_iterator begin_row(size_t row) const {
		assert(_matrix);
		return _matrix->begin_row(row);
	}
	const_row_iterator end_row(size_t row) const {
		assert(_matrix);
		return _matrix->end_row(row);
	}

	const_row_iterator begin_col(size_t col) const {
		return _matrix_transpose.begin_row(col);
	}
	const_row_iterator end_col(size_t col) const {
		return _matrix_transpose.end_row(col);
	}

private:
	ConstSmartPtr<T> _matrix;
	T _matrix_transpose;
};

} // ug

#endif // guard

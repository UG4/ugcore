/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
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


#ifndef __H__UG__COMMON__VARIABLE_ARRAY_H__
#define __H__UG__COMMON__VARIABLE_ARRAY_H__

#include "storage.h"
#include <iostream>
#include <cassert>

namespace ug{

/////////////////////////////////////////////////////////////////////////////////////////////
//	VariableArray1
/**
 * VariableArray1 is a one-dimensional array which supports most of the interface of std::vector
 * (some functions and iterators have to be added).
 * You can use VariableArray1 in DenseVector to get a variable size mathematical Vector.
 * Use storage_traits<type>::is_static to distinguish between fixed and variable vectors.
 * There is a small difference between std::vector and VariableArray1:
 * VariableArray1::capacity == VariableArray1::size, that means there is no capacity buffering.
 * This might be helpful if you have millions of instances of VariableArray1, since you will have
 * a) exactly as many memory allocated as you need
 * b) sizeof(VariableArray1) < sizeof(std::vector) ( 8 vs. 12 on 32-bit systems).
 * \param T type of object in Array (i.e. double)
 */
template<typename T>
class VariableArray1
{
public:
	using value_type = T;
	using size_type = size_t;
	using storage_type = variable_type;

public:
	// 'tors
	VariableArray1();
	explicit VariableArray1(size_type n_);
	VariableArray1(const VariableArray1 &other);

//protected:
	// see Alexandrescu: non-virtual destructors should be protected
	~VariableArray1();

public:
	// capacity
	inline size_type
	size() const;

	inline size_type
	capacity() const;

	inline bool
	resize(size_type n, bool bCopyValues=true);

	inline bool
	reserve(size_type n) const;

	// Element access

	inline const T &
	at (size_type i) const
	{
		// todo: throw if(i >= n)
		return operator [] (i);
	}

	inline T &
	at (size_type i)
	{
		// todo: throw if(i >= n)
		return operator [] (i);
	}

	inline const T &
	operator [] (size_type i) const ;

	inline T &
	operator [] (size_type i) ;

	// output

	template<typename _T>
	friend std::ostream &operator << (std::ostream &out, const VariableArray1<_T> &arr);

protected:
	T *values;
	size_type n;
};

template<typename T>
struct storage_traits1<VariableArray1<T> >
{
	enum {is_static = false};
	enum {static_size = 0};
};

/////////////////////////////////////////////////////////////////////////////////////////////
//	VariableArray2
/**
 * VariableArray2 is a two-dimensional array which supports a similar interface like stl::vector.
 * You can use FixedArray2 in GeMatrix to get a variable size Matrix.
 * Use is_static to distinguish between fixed and variable arrays.
 * \param T type of object in Array (i.e. double)
 * \param T_ordering Ordering of columns/rows. Default is ColMajor. \sa eMatrixOrdering
 */
template<typename T, eMatrixOrdering T_ordering=ColMajor>
class VariableArray2
{
public:
	using value_type = T;
	using size_type = size_t;
	static constexpr eMatrixOrdering ordering = T_ordering;
	enum { is_static=false };
	enum { static_num_rows=0};
	enum { static_num_cols=0};

	using storage_type = variable_type;

public:
	VariableArray2();
	VariableArray2(size_type rows_, size_type cols_);
	VariableArray2(const VariableArray2<T, T_ordering> &other);

//protected:
	// see Alexandrescu: non-virtual destructors should be protected
	~VariableArray2();

public:
	// Capacity
	inline size_type
	num_rows() const;

	inline size_type
	num_cols() const;

	inline bool
	resize(size_type newRows, size_type newCols, bool bCopyValues=true);

	inline size_type
	capacity_num_rows() const { return rows; }

	inline size_type
	capacity_num_cols() const { return cols; };

	inline void
	reserve(size_type nrRows, size_type nrCols) const 	{}

	// Element Access

	inline const T &
	at(size_type r, size_type c) const
	{
		// todo: if(r >= rows || c >= cols) throw
		return operator () (r, c);
	}

	inline T &
	at(size_type r, size_type c)
	{
		// todo: if(r >= rows || c >= cols) throw
		return operator () (r, c);
	}

	inline const T &
	operator () (size_type r, size_type c) const ;

	inline T &
	operator () (size_type r, size_type c) ;

	// output

	template<typename _T, eMatrixOrdering _T_Ordering>
	friend std::ostream &operator << (std::ostream &out, const VariableArray2<_T, _T_Ordering> &arr);

protected:
	T *values;
	size_type rows;
	size_type cols;
};

}

#include "variable_array_impl.h"

#endif
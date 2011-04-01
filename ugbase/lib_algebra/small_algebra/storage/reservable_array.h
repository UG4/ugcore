/*
 *  variable_array.h
 *
 *  Created by Martin Rupp on 21.07.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__UG__COMMON__VARIABLE_ARRAY_H__
#define __H__UG__COMMON__VARIABLE_ARRAY_H__

#include "storage.h"
#include <iostream>
#include <cassert>

namespace ug{

/////////////////////////////////////////////////////////////////////////////////////////////
//	ReservableArray2
/**
 * ReservableArray2 is a two-dimensional array which supports a similar interface like stl::vector.
 * You can use FixedArray2 in GeMatrix to get a variable size Matrix.
 * Use is_static to distinguish between fixed and variable arrays.
 * \param T type of object in Array (i.e. double)
 * \param T_ordering Ordering of columns/rows. Default is ColMajor. \sa eMatrixOrdering
 */
template<typename T, eMatrixOrdering T_ordering=ColMajor>
class ReservableArray2
{
public:
	typedef T value_type;
	typedef size_t size_type;
	static const eMatrixOrdering ordering = T_ordering;
	enum { is_static=false };
	enum { static_num_rows=0};
	enum { static_num_cols=0};
	typedef variable_type storage_type;

public:
	ReservableArray2();
	ReservableArray2(size_type rows_, size_type cols_);
	ReservableArray2(const VariableArray2<T, T_ordering> &other);

//protected:
	// see Alexandrescu: non-virtual destructors should be protected
	~ReservableArray2();

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

	inline bool
	reserve(size_type nrRows, size_type nrCols) const 	{ return; }

	// Element Access

	inline const T &
	at(size_type r, size_type c) const
	{
		// todo: if(r >= rows || c >= cols) throw
		return operator()(r, c);
	}

	inline T &
	at(size_type r, size_type c)
	{
		// todo: if(r >= rows || c >= cols) throw
		return operator()(r, c);
	}

	inline const T &
	operator()(size_type r, size_type c) const ;

	inline T &
	operator()(size_type r, size_type c) ;

	// output

	template<typename _T, eMatrixOrdering _T_Ordering>
	friend std::ostream &operator << (std::ostream &out, const ReservableArray2<_T, _T_Ordering> &arr);

protected:
	T *values;
	size_type rows;
	size_type cols;
	size_type arraySize;
};

}

#include "variable_array_impl.h"

#endif // #define __H__UG__COMMON__VARIABLE_ARRAY_H__

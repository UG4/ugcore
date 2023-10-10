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


#ifndef __H__UG__COMMON__FIXED_ARRAY_IMPL_H__
#define __H__UG__COMMON__FIXED_ARRAY_IMPL_H__

#include "storage.h"
#include "fixed_array.h"

namespace ug{
// 'tors

template<typename T, size_t n>
FixedArray1<T, n>::FixedArray1()
{
}


template<typename T, size_t nT>
FixedArray1<T, nT>::FixedArray1(size_t n)
{
	assert(n == nT);
}

template<typename T, size_t n>
FixedArray1<T, n>::FixedArray1(const FixedArray1<T, n> &other)
{
	for(size_type i=0; i<n; i++)
		values[i] = other[i];
}

template<typename T, size_t n>
FixedArray1<T, n>::~FixedArray1()
{
}

// capacity

template<typename T, size_t n>
inline size_t
FixedArray1<T, n>::size() const
{
	return n;
}


template<typename T, size_t n>
inline bool
FixedArray1<T, n>::resize(size_t newN, bool bCopyValues)
{
	assert(newN == n);
	return newN == n;
}


template<typename T, size_t n>
inline bool
FixedArray1<T, n>::reserve(size_t newN) const
{
	assert(newN <= n);
	return newN <= n;
}


// Element access
template<typename T, size_t n>
inline T &
FixedArray1<T, n>::operator[](size_t i)
{
	assert(i<n);
	return values[i];
}

template<typename T, size_t n>
inline const T &
FixedArray1<T, n>::operator[](size_t i) const
{
	assert(i<n);
	return values[i];
}

// output

template<typename T, size_t n>
std::ostream &operator << (std::ostream &out, const FixedArray1<T, n> &arr)
{
	out << "FixedArray (n=" << n << ") [ ";
	for(typename FixedArray1<T, n>::size_type i=0; i<arr.size(); i++)
		out << arr[i] << " ";
	out << "]";
	return out;
}


////////////////////////////////////////////////////////////////////////////////

// 'tors

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
FixedArray2<T, rowsT, colsT, T_ordering>::FixedArray2()
{
}


template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
FixedArray2<T, rowsT, colsT, T_ordering>::FixedArray2(size_t rows, size_t cols)
{
	assert(rows == rowsT && cols == colsT);
}

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
FixedArray2<T, rowsT, colsT, T_ordering>::FixedArray2(const FixedArray2<T, rowsT, colsT, T_ordering> &other)
{
	for(size_type i=0; i<rowsT*colsT; i++)
		values[i] = other.values[i];
}

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
FixedArray2<T, rowsT, colsT, T_ordering>::~FixedArray2()
{
}

// Capacity

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
inline size_t
FixedArray2<T, rowsT, colsT, T_ordering>::num_rows() const
{
	return rowsT;
}

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
inline size_t
FixedArray2<T, rowsT, colsT, T_ordering>::num_cols() const
{
	return colsT;
}


template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
inline bool
FixedArray2<T, rowsT, colsT, T_ordering>::resize(size_t newRows, size_t newCols, bool bCopyValues)
{
	assert(newCols == colsT && newRows == rowsT);
	return newCols == colsT && newRows == rowsT;
}



// Element access

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
inline T &
FixedArray2<T, rowsT, colsT, T_ordering>::operator()(size_t r, size_t c)
{
	assert(r<rowsT);
	assert(c<colsT);
	if(T_ordering == RowMajor)
		return values[c+r*colsT];
	else
		return values[r+c*rowsT];
}

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
inline const T &
FixedArray2<T, rowsT, colsT, T_ordering>::operator()(size_t r, size_t c) const
{
	assert(r<rowsT);
	assert(c<colsT);
	if(T_ordering == RowMajor)
		return values[c+r*colsT];
	else
		return values[r+c*rowsT];
}

/*
template<typename T, size_t rowsT, size_t colsT>
T &
FixedArray2<T, rowsT, colsT, ColMajor>::operator()(size_t r, size_t c)
{
	assert(r>=0 && r<rowsT);
	assert(c>=0 && c<colsT);
	return values[c+r*colsT];
}

template<typename T, size_t rowsT, size_t colsT>
T &
FixedArray2<T, rowsT, colsT, RowMajor>::operator()(size_t r, size_t c)
{
	assert(r>=0 && r<rowsT);
	assert(c>=0 && c<colsT);
	return values[r+c*rowsT];
}


template<typename T, size_t rowsT, size_t colsT>
const T &
FixedArray2<T, rowsT, colsT, ColMajor>::operator()(size_t r, size_t c) const
{
	assert(r>=0 && r<rowsT);
	assert(c>=0 && c<colsT);
	return values[c+r*colsT];
}

template<typename T, size_t rowsT, size_t colsT>
const T &
FixedArray2<T, rowsT, colsT, RowMajor>::operator()(size_t r, size_t c) const
{
	assert(r>=0 && r<rowsT);
	assert(c>=0 && c<colsT);
	return values[r+c*rowsTT];
}*/


// output

template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering>
std::ostream &operator << (std::ostream &out, const FixedArray2<T, rowsT, colsT, T_ordering> &arr)
{
	//out << "FixedArray2 (" << rowsT << "x" << colsT << "), " << ((T_ordering == ColMajor) ? "ColMajor" : "RowMajor") << endl;
	out << "[";
	typedef typename FixedArray2<T, rowsT, colsT, T_ordering>::size_type size_type;
	for(size_type r=0; r<arr.num_rows(); r++)
	{
		for(size_type c=0; c<arr.num_cols(); c++)
			out << arr(r, c) << " ";
		if(r != arr.num_rows()-1) out << "| ";
	}
	out << "]";

	return out;
}

}
#endif  // __H__UG__COMMON__FIXED_ARRAY_IMPL_H__

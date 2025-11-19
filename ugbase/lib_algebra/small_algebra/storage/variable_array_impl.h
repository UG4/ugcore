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


#ifndef __H__UG__COMMON__VARIABLE_ARRAY_IMPL_H__
#define __H__UG__COMMON__VARIABLE_ARRAY_IMPL_H__

#include "storage.h"
#include "variable_array.h"
#include "common/common.h"
#include <algorithm> // for min
#include <cstring>

namespace ug{

// 'tors
template<typename T>
VariableArray1<T>::VariableArray1()
{
	n = 0;
	values = nullptr;
}


template<typename T>
VariableArray1<T>::VariableArray1(size_t n_)
{
	n = 0;
	values = nullptr;
	resize(n_, false);
}

template<typename T>
VariableArray1<T>::VariableArray1(const VariableArray1<T> &other)
{
	if(this == &other) return;
	n = 0;
	values = nullptr;
	resize(other.size(), false);
	for(size_type i=0; i<n; i++)
		values[i] = other[i];
}

template<typename T>
VariableArray1<T>::~VariableArray1()
{
	if(values) { delete[] values; values = nullptr; }
	n = 0;
}


// Capacity

template<typename T>
inline size_t
VariableArray1<T>::size() const
{
	return n;
}

template<typename T>
bool
VariableArray1<T>::resize(size_t newN, bool bCopyValues)
{
	if(newN == n) return true;

	if(newN <= 0)
	{
		if(values) delete[] values;
		values = nullptr;
		n = 0;
		return true;
	}
	value_type *new_values = new T[newN];
	UG_ASSERT(new_values != nullptr, "out of memory");
	if(new_values == nullptr) return false;
	memset(reinterpret_cast<void *> (new_values), 0, sizeof(T)*newN); // todo: think about that

	if(bCopyValues)
	{
		/*
		if(storage_traits<value_type>::is_static)
		{
			for(int i=0; i<n; i++)
				new_values[i] = values[i];
		}
		else {
		// we are using swap to avoid re-allocations
		 */
		size_t minN = std::min(n, newN);
		for(size_t i=0; i<minN; i++)
			std::swap(new_values[i], values[i]);
	}

	if(values) delete[] values;
	values = new_values;
	n = newN;
	return true;
}

template<typename T>
inline size_t
VariableArray1<T>::capacity() const
{
	return n;
}


// use stl::vector if you want to use reserve
template<typename T>
inline bool
VariableArray1<T>::reserve(size_t newCapacity) const
{
	return true;
}


// Element access

template<typename T>
T &
VariableArray1<T>::operator[](size_t i)
{
	assert(values);
	assert(i<n);
	return values[i];
}

template<typename T>
const T &
VariableArray1<T>::operator[](size_t i) const
{
	assert(values);
	assert(i<n);
	return values[i];
}

template<typename T>
std::ostream &operator << (std::ostream &out, const VariableArray1<T> &arr)
{
	//out << "VariableArray (n=" << arr.size() << ") [ ";
	for(size_t i=0; i<arr.size(); i++)
		out << arr[i] << " ";
	out << "]";
	return out;
}


////////////////////////////////////////////////////////////////////////////////


template<typename T, eMatrixOrdering T_ordering>
VariableArray2<T, T_ordering>::VariableArray2()
{
	values = nullptr;
	rows = 0;
	cols = 0;
}


template<typename T, eMatrixOrdering T_ordering>
VariableArray2<T, T_ordering>::VariableArray2(size_t rows, size_t cols)
{
	values = nullptr;
	rows = 0;
	cols = 0;
	resize(rows, cols);
}

template<typename T, eMatrixOrdering T_ordering>
VariableArray2<T, T_ordering>::VariableArray2(const VariableArray2<T, T_ordering> &other)
{
	if(this == &other) return;
	values = nullptr;
	rows = 0;
	cols = 0;
	resize(other.num_rows(), other.num_cols(), false);
	for(size_type i=0; i<rows*cols; i++)
		values[i] = other.values[i];
}

template<typename T, eMatrixOrdering T_ordering>
VariableArray2<T, T_ordering>::~VariableArray2()
{
	if(values) { delete[] values; values = nullptr; }
	rows = cols = 0;
}

// Capacity

template<typename T, eMatrixOrdering T_ordering>
size_t
VariableArray2<T, T_ordering>::num_rows() const
{
	return rows;
}

template<typename T, eMatrixOrdering T_ordering>
size_t
VariableArray2<T, T_ordering>::num_cols() const
{
	return cols;
}


template<typename T, eMatrixOrdering T_ordering>
bool
VariableArray2<T, T_ordering>::resize(size_t newRows, size_t newCols, bool bCopyValues)
{
	if(newRows == rows && newCols == cols) return true;

	if(newRows == 0 && newCols == 0)
	{
		rows = cols = 0;
		if(values) delete[] values;
		values = nullptr;
		return true;
	}

	value_type *new_values = new T[newRows*newCols];
	memset(reinterpret_cast<void *> (new_values), 0, sizeof(T)*newRows*newCols); // todo: think about that
	UG_ASSERT(new_values != nullptr, "out of memory");
	if(new_values==nullptr) return false;
	/*
	if(storage_traits<value_type>::is_static)
	{
		...
	}
	else {

	 */
	if(bCopyValues)
	{
		size_t minRows = std::min(rows, newRows);
		size_t minCols = std::min(cols, newCols);

		// we are using swap to avoid re-allocations
		if(T_ordering==RowMajor)
			for(size_t r=0; r<minRows; r++)
				for(size_t c=0; c<minCols; c++)
					std::swap(new_values[c+r*newCols], values[c+r*cols]);
		else
			for(size_t r=0; r<minRows; r++)
				for(size_t c=0; c<minCols; c++)
					std::swap(new_values[r+c*newRows], values[r+c*rows]);
	}

	if(values) delete[] values;
	rows = newRows;
	cols = newCols;
	values = new_values;
	return true;
}


template<typename T, eMatrixOrdering T_ordering>
T &
VariableArray2<T, T_ordering>::operator()(size_t r, size_t c)
{
	UG_ASSERT(r<rows, "r = " << r << ", rows = " << rows);
	UG_ASSERT(c<cols, "c = " << c << ", cols = " << cols);
	if(T_ordering==RowMajor)
		return values[c+r*cols];
	else
		return values[r+c*rows];
}

template<typename T, eMatrixOrdering T_ordering>
const T &
VariableArray2<T, T_ordering>::operator()(size_t r, size_t c) const
{
	UG_ASSERT(r<rows, "r = " << r << ", rows = " << rows);
	UG_ASSERT(c<cols, "c = " << c << ", cols = " << cols);
	if(T_ordering==RowMajor)
		return values[c+r*cols];
	else
		return values[r+c*rows];
}

template<typename T, eMatrixOrdering T_ordering>
std::ostream &operator << (std::ostream &out, const VariableArray2<T, T_ordering> &arr)
{
	out << "[ ";
	//out << "VariableArray2 (" << arr.num_rows() << "x" << arr.num_cols() << "), " << ((T_ordering == ColMajor) ? "ColMajor" : "RowMajor") << endl;
	using size_type  = size_t ;
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
#endif

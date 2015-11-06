/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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
#include <algorithm> // for min
#include <cstring>

namespace ug{


////////////////////////////////////////////////////////////////////////////////


template<typename T, eMatrixOrdering T_ordering>
ReservableArray2<T, T_ordering>::ReservableArray2()
{
	values = NULL;
	rows = 0;
	cols = 0;
	arraySize = 0;
}


template<typename T, eMatrixOrdering T_ordering>
ReservableArray2<T, T_ordering>::ReservableArray2(size_t rows, size_t cols)
{
	values = NULL;
	rows = 0;
	cols = 0;
	arraySize = 0;
	resize(rows, cols);
}

template<typename T, eMatrixOrdering T_ordering>
ReservableArray2<T, T_ordering>::ReservableArray2(const ReservableArray2<T, T_ordering> &other)
{
	if(this == &other) return;
	if(values) { delete[] values; values = NULL; }
	rows = 0;
	cols = 0;
	arraySize = 0;
	resize(other.num_rows(), other.num_cols());
	for(size_type i=0; i<rows*cols; i++)
		values[i] = other.values[i];
}

template<typename T, eMatrixOrdering T_ordering>
ReservableArray2<T, T_ordering>::~ReservableArray2()
{
	if(values) { delete[] values; values = NULL; }
	rows = cols = 0;
}

// Capacity

template<typename T, eMatrixOrdering T_ordering>
size_t
ReservableArray2<T, T_ordering>::num_rows() const
{
	return rows;
}

template<typename T, eMatrixOrdering T_ordering>
size_t
ReservableArray2<T, T_ordering>::num_cols() const
{
	return cols;
}


template<typename T, eMatrixOrdering T_ordering>
bool
ReservableArray2<T, T_ordering>::resize(size_t newRows, size_t newCols, bool bCopyValues)
{
	assert(newRows >= 0 && newCols >= 0);
	if(newRows == rows && newCols == cols) return true;

	if(newRows == 0 && newCols == 0)
	{
		rows = cols = 0;
		if(values) delete[] values;
		values = NULL;
		return true;
	}

	if(bCopyValues)
	{
		size_t minRows = std::min(rows, newRows);
		size_t minCols = std::min(cols, newCols);
		if(newRows*newCols > arraySize)
		{
			value_type *new_values = new T[newRows*newCols];
			arraySize = newRows*newCols;
			memset(new_values, 0, sizeof(T)*newRows*newCols); // todo: think about that
			UG_ASSERT(new_values != NULL, "out of memory");
			if(new_values==NULL) return false;
			/*
			if(storage_traits<value_type>::is_static)
			{
				...
			}
			else {

			 */

			// we are using swap to avoid re-allocations
			if(T_ordering==RowMajor)
				for(size_t r=0; r<minRows; r++)
					for(size_t c=0; c<minCols; c++)
						std::swap(new_values[c+r*newCols], values[c+r*cols]);
			else
				for(size_t r=0; r<minRows; r++)
					for(size_t c=0; c<minCols; c++)
						std::swap(new_values[r+c*newRows], values[r+c*rows]);

			if(values) delete[] values;
		}
		else
		{
			if(T_ordering==RowMajor)
			{
				if(newCols < cols)
				{
					for(size_t r=1; r<minRows; r++)
						for(size_t c=0; c<minCols; c++)
							std::swap(values[c+r*newCols], values[c+r*cols]);

				}
				else if(newCols > cols)
				{
					for(size_t r=minRows-1; r>0; r--)
					{
						size_t c=newCols-1;
						do
							std::swap(values[c+r*newCols], values[c+r*cols]);
						while(--c > 0);
					}
				}
				else
				{
					// nothing to do
				}
			}
			else
			{
				if(newRows < rows)
				{
					for(size_t c=1; c<minCols; c++)
						for(size_t r=0; r<minRows; r++)
							std::swap(values[r+c*newRows], values[r+c*rows]);

				}
				else if(newRows > rows)
				{
					for(size_t c=minCols-1; c>0; c--)
					{
						size_t r=newRows-1;
						do
							std::swap(values[r+c*newRows], values[r+c*rows]);
						while(--r > 0);
					}
				}
				else
				{
					// nothing to do
				}
			}

		}
	}
	rows = newRows;
	cols = newCols;
	values = new_values;
	return true;
}


template<typename T, eMatrixOrdering T_ordering>
T &
ReservableArray2<T, T_ordering>::operator()(size_t r, size_t c)
{
	assert(r>=0 && r<rows);
	assert(c>=0 && c<cols);
	if(T_ordering==RowMajor)
		return values[c+r*cols];
	else
		return values[r+c*rows];
}

template<typename T, eMatrixOrdering T_ordering>
const T &
ReservableArray2<T, T_ordering>::operator()(size_t r, size_t c) const
{
	assert(r>=0 && r<rows);
	assert(c>=0 && c<cols);
	if(T_ordering==RowMajor)
		return values[c+r*cols];
	else
		return values[r+c*rows];
}

template<typename T, eMatrixOrdering T_ordering>
std::ostream &operator << (std::ostream &out, const ReservableArray2<T, T_ordering> &arr)
{
	out << "[ ";
	//out << "ReservableArray2 (" << arr.num_rows() << "x" << arr.num_cols() << "), " << ((T_ordering == ColMajor) ? "ColMajor" : "RowMajor") << endl;
	typedef size_t size_type;
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
#endif // __H__UG__COMMON__VARIABLE_ARRAY_IMPL_H__

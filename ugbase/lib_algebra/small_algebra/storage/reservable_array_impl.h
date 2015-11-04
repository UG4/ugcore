
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

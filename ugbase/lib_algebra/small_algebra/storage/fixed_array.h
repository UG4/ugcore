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


#ifndef __H__UG__COMMON__FIXED_ARRAY_H__
#define __H__UG__COMMON__FIXED_ARRAY_H__

#include "storage.h"
#include <iostream>
#include <cassert>

/////////////////////////////////////////////////////////////////////////////////////////////
//	FixedArray1
namespace ug{

/**
 * FixedArray1 is a one-dimensional array which supports most of the interface of std::vector.
 * (some functions and iterators have to be added).
 * You can use FixedArray1 in DenseVector to get a fixed size mathematical Vector.
 * Use storage_traits<type>::is_static to distinguish between fixed and variable vectors. Functions like
 * resize(newN) assert N == newN. So FixedArray1<T, N>::size is always N.
 * \tparam T type of object in Array (i.e. double)
 * \tparam n number of elements in the fixed array (T values[n]; )
 */
template<typename T, size_t n>
class FixedArray1
{
public:
	using value_type = T;
	using size_type = size_t;
	using storage_type = static_type;

public:
	FixedArray1() = default;
	explicit FixedArray1(size_type n_);
	FixedArray1(const FixedArray1 &other);

//protected:
	// see Alexandrescu: non-virtual destructors should be protected
	~FixedArray1() = default;

public:
	// capacity
	inline size_type
	size() const;

	inline size_type
	capacity() const;

	inline bool
	resize(size_type newN, bool bCopyValues=true);

	inline bool
	reserve(size_type newN) const;

	// Element access
	inline const T &
	at(size_type i) const
	{
		return operator [] (i);
	}

	inline T &
	at(size_type i)
	{
		return operator [] (i);
	}

	inline const T &
	operator [] (size_type i) const ;

	inline T &
	operator [] (size_type i) ;

	// output
	template<typename _T, size_type _n>
	friend
	std::ostream &
	operator << (std::ostream &out, const FixedArray1<_T, _n> &arr);

protected:
	T values[n];
};

template<typename T, size_t N>
struct storage_traits1<FixedArray1<T, N> >
{
	enum {is_static = true};
	enum {static_size = N};
};


/////////////////////////////////////////////////////////////////////////////////////////////
//	FixedArray2
/**
 * FixedArray2 is a two-dimensional array which supports a similar interface like stl::vector.
 * You can use FixedArray2 in GeMatrix to get a fixed size Matrix.
 * Use is_static to distinguish between fixed and variable arrays.
 * \tparam T type of object in Array (i.e. double)
 * \tparam colsT fixed number of columns
 * \tparam rowsT fixed number of rows
 * \tparam T_ordering Ordering of columns/rows. Default is ColMajor. \sa eMatrixOrdering
 */
template<typename T, size_t rowsT, size_t colsT, eMatrixOrdering T_ordering=ColMajor>
class FixedArray2
{
public:
	using value_type = T;
	using size_type = size_t;
	static constexpr eMatrixOrdering ordering = T_ordering;
	enum { is_static=true};
	enum { static_num_rows=rowsT};
	enum { static_num_cols=colsT};

	using storage_type = static_type;

public:
	FixedArray2() = default;
	FixedArray2(size_type rows, size_type cols);
	FixedArray2(const FixedArray2 &other);

//protected:
	// see Alexandrescu: non-virtual destructors should be protected
	~FixedArray2() = default;

public:
	// capacity
	inline size_type
	num_rows() const;

	inline size_type
	num_cols() const;

	inline bool
	resize(size_type newRows, size_type newCols, bool bCopyValues=true);

	inline size_type
	capacity_num_rows() const { return rowsT; }

	inline size_type
	capacity_num_cols() const { return colsT; };

	inline bool
	empty() const;

	inline bool
	reserve(size_type nrRows, size_type nrCols) const
	{
		assert(nrRows == rowsT && nrCols == colsT);
		return nrRows == rowsT && nrCols == colsT;
	}


	// Element access
	inline const T &
	at(size_type r, size_type c) const
	{
		// todo: if(r >= rowsT || c >= colsT) throw
		return operator () (r, c);
	}

	inline T &
	at(size_type r, size_type c)
	{
		// todo: if(r >= rowsT || c >= colsT) throw
		return operator () (r, c);
	}

	inline const T &
	operator () (size_type r, size_type c) const ;

	inline T &
	operator () (size_type r, size_type c) ;

	// output
	template<typename a, size_type b, size_type c, eMatrixOrdering d>
	friend
	std::ostream &
	operator << (std::ostream &out, const FixedArray2<a, b, c, d> &arr);

protected:
	T values[rowsT*colsT];
};


}
#include "fixed_array_impl.h"
// #include "fixed_array_specialization.h"

#endif

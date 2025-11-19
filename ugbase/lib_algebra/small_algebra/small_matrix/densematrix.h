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


#ifndef __H__UG__COMMON__DENSEMATRIX_H__
#define __H__UG__COMMON__DENSEMATRIX_H__

#include <iostream>
#include <cassert>
#include "../storage/storage.h"
#include "densevector.h"

namespace ug{

/// \addtogroup small_algebra
/// \{

/////////////////////////////////////////////////////////////////////////////////////////////
//	DenseMatrix
/**
 * DenseMatrix is a mathematical matrix class, which inherits its storage behaviour
 * (fixed/variable size, RowMajor/ColMajor ordering) from TStorage.
 * \param TStorage storage policy with interface of VariableArray2.
 * \sa FixedArray2, VariableArray2, eMatrixOrdering
 */
template<typename TStorage>
class DenseMatrix : public TStorage
{
public:
	using value_type = typename TStorage::value_type;
	using size_type = typename TStorage::size_type;
	static constexpr eMatrixOrdering ordering = TStorage::ordering;
	enum { is_static = TStorage::is_static};
	enum { static_num_rows = TStorage::static_num_rows};
	enum { static_num_cols = TStorage::static_num_cols};

	using this_type = DenseMatrix<TStorage>;
	using base = TStorage;
	using base::operator ();
	using base::at;
	using base::num_rows;
	using base::num_cols;
	using base::resize;
	
public:
	// 'tors
	DenseMatrix();
	DenseMatrix(double val);
	DenseMatrix(const this_type &rhs);

	//~DenseMatrix() {} // dont implement a destructor, since ~base may not be virtual

public:
	// matrix assignment operators

	
////// =
	inline
	this_type &
	operator = (const this_type &t);

	template<typename T>
	inline
	this_type &
	operator = (const T &t);
	
	inline
	this_type &
	operator = (double alpha);
	
////// +=
	template<typename T>
	inline	
	this_type &
	operator += (const T &t);	
	
	inline
	this_type &
	operator += (double alpha);	
		
////// -=
	template<typename T>
	inline	
	this_type &
	operator-=(const T &t);
	
	inline
	this_type &
	operator-=(double alpha);

	
////// *=
	inline
	this_type&
	operator*=(double alpha);
	
	inline
	this_type&
	operator*=(const this_type &mat);

////// /=
	inline
	this_type&
	operator/=(double alpha);
	
	inline
	this_type&
	operator /= (this_type &other);
	
	
////// +
	inline
	this_type
	operator+(const this_type &other ) const;
	
////// -
	inline
	this_type
	operator-(const this_type &other ) const;
	
	inline
	this_type
	T()  const;

////// unary -
	inline
	this_type
	operator-() const;

////// *

	
	// multiply
	template<typename TStorage2>
	DenseVector<TStorage2>
	operator * (const DenseVector<TStorage2> &vec) const;

	inline
	this_type
	operator * (const this_type &mat) const;

	inline
	this_type
	operator * (double alpha ) const;

///// /
	this_type
	operator / (this_type &other);

// compare operators
////// ==
	inline
	bool
	operator == (double t) const;

	
	template<typename TStorage2>
	inline
	bool
	operator == (const DenseMatrix<TStorage2> &t) const;
	
///// !=
	template<typename TStorage2>
	inline
	bool
	operator != (const TStorage2 &t) const;

////// other
	inline
	const value_type &
	entry(size_type r, size_type c) const
	{
		return operator()(r,c);
	}

	inline
	value_type &
	entry(size_type r, size_type c)
	{
		return operator()(r,c);
	}

	
	template<typename TStorage2>
	void
	subassign(size_t r, size_t c, const TStorage2 &t)
	{
		UG_ASSERT(r+t.num_rows() <= num_rows() && c+t.num_cols() <= num_cols(), "");
		for(size_t r1=0; r1<t.num_rows(); r1++)
			for(size_t c1=0; c1<t.num_cols(); c1++)
				entry(r+r1, c+c1) = t(r1, c1);			
	}
	
	void
	subassign(size_t r, size_t c, const value_type &t)
	{
		entry(r, c) = t;
	}

	void maple_print(const char *name);
};

template<typename TStorage>
DenseMatrix<TStorage> operator *(number a, const DenseMatrix<TStorage> &b)
{
	return b*a;
}


template<typename TStorage>
DenseMatrix<TStorage> MatrixTranspose(const DenseMatrix<TStorage> &A)
{
	DenseMatrix<TStorage> At;
	At.resize(A.num_cols(), A.num_rows());
	for(size_t r=0; r < A.num_rows(); r++)
		for(size_t c=0; c < A.num_cols(); c++)
			At(c,r) = A(r, c);
	return At;
}


inline double MatrixTranspose(const double &b)
{return b;}
// end group small_algebra
/// \}

}

#include "densematrix_impl.h"
#include "densematrix_operations.h"

#endif
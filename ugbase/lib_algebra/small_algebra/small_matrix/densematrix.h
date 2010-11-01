/**
 * \file densematrix.h
 *
 * \author Martin Rupp
 *
 * \date 21.07.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */


#ifndef __H__UG__COMMON__DENSEMATRIX_H__
#define __H__UG__COMMON__DENSEMATRIX_H__

#include <iostream>
#include <cassert>
#include "../storage/storage.h"

namespace ug{

/////////////////////////////////////////////////////////////////////////////////////////////
//	DenseMatrix
/**
 * DenseMatrix is a mathematical matrix class, which inheritates its storage behaviour
 * (fixed/variable size, RowMajor/ColMajor ordering) from TStorage.
 * \param TStorage storage policy with interface of VariableArray2.
 * \sa FixedArray2, VariableArray2, eMatrixOrdering
 */
template<typename TStorage>
class DenseMatrix : public TStorage
{
public:
	typedef typename TStorage::value_type value_type;
	typedef typename TStorage::size_type size_type;
	static const eMatrixOrdering ordering = TStorage::ordering;
	enum { is_static = TStorage::is_static};
	enum { static_num_rows = TStorage::static_num_rows};
	enum { static_num_cols = TStorage::static_num_cols};

	typedef DenseMatrix<TStorage> this_type;
	typedef TStorage base;
	using base::operator ();
	using base::at;
	using base::num_rows;
	using base::num_cols;

public:
	// 'tors
	DenseMatrix();
	DenseMatrix(const this_type &rhs);

	//~DenseMatrix() {} // dont implement a destructor, since ~base may not be virtual

public:
	// matrix assignment operators
	inline this_type &
	operator =  (const this_type &rhs);

	inline this_type &
	operator += (const this_type &rhs);

	inline this_type &
	operator -= (const this_type &rhs);

	// alpha operators
	template<typename T>
	inline this_type &
	operator=(const T &alpha);

	inline this_type &
	operator+=(const value_type &alpha);

	inline this_type &
	operator-=(const value_type &alpha);

	template<typename T>
	inline this_type &
	operator*=(const T &alpha);

	inline this_type &
	operator/=(const value_type &alpha);

	// compare operators
	template<typename T>
	inline bool
	operator == (const T &t) const;

	template<typename T>
	inline bool
	operator == (const DenseMatrix<T> &t) const;

	template<typename T>
	inline bool
	operator != (const T &t) const;


	inline const value_type &
	entry(size_type r, size_type c) const
	{
		return operator()(r,c);
	}

	inline value_type &
	entry(size_type r, size_type c)
	{
		return operator()(r,c);
	}
};

}

#include "densematrix_impl.h"
#include "densematrix_operations.h"


#endif // __H__UG__COMMON__DENSEMATRIX_H__

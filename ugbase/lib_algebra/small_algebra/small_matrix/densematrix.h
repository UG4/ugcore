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

	////////////////////////////////////////////////////////////////////
	// this will be removed soon


	this_type operator + (const this_type &other ) const
	{
		UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), "");
		this_type erg;
		erg.resize(num_rows(), num_cols());
		for(size_t r=0; r<num_rows(); r++)
			for(size_t c=0; c<num_cols(); c++)
				erg(r, c) = entry(r, c) + other(r,c);
		return erg;
	}

	this_type operator - (const this_type &other ) const
	{
		UG_ASSERT(num_rows() == other.num_rows() && num_cols() == other.num_cols(), "");
		this_type erg;
		erg.resize(num_rows(), num_cols());

		for(size_t r=0; r<num_rows(); r++)
			for(size_t c=0; c<num_cols(); c++)
				erg(r, c) = entry(r, c) - other(r,c);
		return erg;
	}

// multiply
	this_type operator * (const this_type &other ) const
	{
		// that aint 100% correct
		UG_ASSERT(num_cols() == other.num_rows(), "");

		this_type erg;
		erg.resize(num_rows(), other.num_cols());

		for(size_t r=0; r < num_rows(); r++)
			for(size_t c=0; c < other.num_cols(); c++)
			{
				for(size_t i=0; i < num_cols(); i++)
					AddMult(erg(r,c), at(r, i), other.at(i, c));
			}
		return erg;
	}



	this_type &operator /= (this_type &other)
	{
		this_type tmp = other;
		Invert(tmp);

		(*this) = (*this) * tmp;
		return *this;
	}

	this_type operator / (this_type &other)
	{
		this_type tmp = other;
		Invert(tmp);

		return (*this) * tmp;
	}
};

template<typename TStorage>
DenseMatrix<TStorage> operator *(number a, const DenseMatrix<TStorage> &b)
{
	return b*a;
}

}

#include "densematrix_impl.h"
#include "densematrix_operations.h"


#endif // __H__UG__COMMON__DENSEMATRIX_H__

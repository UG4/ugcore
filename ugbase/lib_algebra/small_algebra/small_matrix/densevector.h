/**
 * \file densevector.h
 *
 * \author Martin Rupp
 *
 * \date 21.07.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__COMMON__DENSEVECTOR_H__
#define __H__UG__COMMON__DENSEVECTOR_H__

#include <iostream>
#include <cassert>
#include "../storage/storage.h"

namespace ug{


/////////////////////////////////////////////////////////////////////////////////////////////
//	DenseVector
/**
 * DenseVector is a one-dimensional mathematical vector of objects TStorage::value_type.
 * Use DenseVector with FixedArray1, VariableArray1 or stl::vector. Depending on
 * TStorage, DenseVector is of fixed size or variable size, and inheritates the interface
 * or TStorage.
 * \param TStorage storage policy with interface of VariableArray1.
 * \sa FixedArray1, VariableArray1
 */
template<typename TStorage>
class DenseVector : public TStorage
{
public:
	typedef typename TStorage::value_type value_type;
	typedef typename TStorage::size_type size_type;

	// use traits so we are able to use std::vector
	enum { is_static = storage_traits1<TStorage>::is_static};
	enum { static_size = storage_traits1<TStorage>::static_size};

	typedef DenseVector<TStorage> this_type;

	// use the interface of TStorage
	typedef TStorage base;
	using base::operator [];
	using base::at;
	using base::size;
	using base::resize;

public:
	// constructors
	DenseVector();
	DenseVector(const DenseVector<TStorage> &rhs);

	// ~DenseVector(); // dont implement a destructor, since ~base may not be virtual

	// operations with vectors
	inline this_type &
	operator= (const this_type &rhs);


	inline this_type &
	operator+= (const this_type &rhs);

	inline this_type &
	operator-= (const this_type &rhs);

	// operations with scalars
	template<typename T>
	inline this_type&
	operator=  (const T& alpha);

	inline this_type&
	operator+= (const value_type& alpha);

	inline this_type&
	operator-= (const value_type& alpha);

	template<typename T>
	inline this_type&
	operator *= (const T &alpha);

	inline this_type&
	operator/= (const value_type& alpha);


	value_type &entry(size_t i) const
	{
		return operator[] (i);
	}
	value_type &entry(size_t i)
	{
		return operator[] (i);
	}

	template<typename Type>
	DenseVector<TStorage> &
	assign(const Type &t);

	void maple_print(const char *name);

};

}
#include "densevector_impl.h"

#endif // __H__UG__COMMON__DENSEVECTOR_H__

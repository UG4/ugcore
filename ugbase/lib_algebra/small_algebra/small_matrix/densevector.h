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

/// \addtogroup small_algebra
/// \{

template<typename T>
class TE_TRANSPOSED
{
public:
	typedef typename T::value_type value_type;
	TE_TRANSPOSED(const T&_t) : t(_t) {}
	inline size_t num_rows() const
	{
		return t.num_cols();
	}
	
	inline size_t num_cols() const
	{
		return t.num_rows();
	}
	
	const value_type &operator() (size_t r, size_t c) const
	{
		return t(c, r);
	}
	
	value_type &operator() (size_t r, size_t c)
	{
		return t(c, r);
	}	
private:	
	const T &t;
};

template<typename T>
inline TE_TRANSPOSED<T> te_transpose(const T &t)
{
	return TE_TRANSPOSED<T>(t);
}

inline const double &te_transpose(const double &t)
{
	return t;
}

inline double &te_transpose(double &t)
{
	return t;
}
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
	DenseVector(double alpha=0.0);
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
	
	template<typename T>
	inline this_type&
	operator=  (const double & alpha)
	{
		for(size_t i=0; i<size(); i++)
			entry(i) = alpha;
		return *this;
	}

	inline this_type&
	operator+= (const value_type& alpha);

	inline this_type&
	operator-= (const value_type& alpha);

	template<typename T>
	inline this_type&
	operator *= (const T &alpha);

	inline this_type&
	operator/= (const value_type& alpha);


	bool operator > (double d) const
	{
		for(size_t i=0; i<size(); i++)
			if(entry(i) <= d) return false;
		return true;
	}

	bool operator < (double d) const
	{
		for(size_t i=0; i<size(); i++)
			if(entry(i) >= d) return false;
		return true;
	}

	bool operator >= (double d) const
	{
		for(size_t i=0; i<size(); i++)
			if(entry(i) < d) return false;
		return true;
	}

	bool operator <= (double d) const
	{
		for(size_t i=0; i<size(); i++)
			if(entry(i) > d) return false;
		return true;
	}
	bool operator == (double d) const
	{
		for(size_t i=0; i<size(); i++)
			if(entry(i) != d) return false;
		return true;
	}

	bool operator != (double d) const
	{
		for(size_t i=0; i<size(); i++)
			if(entry(i) != d) return true;
		return false;
	}
	////////////////////////////////////////////////////////////////////
	// this will be removed soon


	this_type operator + (const this_type &other ) const
	{
		UG_ASSERT(size() == other.size(), "");
		this_type erg;
		erg.resize(size());
		for(size_t i=0; i<size(); i++)
			erg[i] = entry(i) + other[i];
		return erg;
	}

	this_type operator - (const this_type &other ) const
	{
		UG_ASSERT(size() == other.size(), "");
		this_type erg;
		erg.resize(size());
		for(size_t i=0; i<size(); i++)
			erg[i] = entry(i) - other[i];
		return erg;
	}

// multiply
	this_type operator * (double alpha ) const
	{
		this_type erg;
		erg.resize(size());
		for(size_t i=0; i<size(); i++)
			erg[i] = alpha*entry(i);
		return erg;
	}
	
	this_type operator - () const
	{
		this_type erg;
		erg.resize(size());
		for(size_t i=0; i<size(); i++)
			erg[i] *= -1.0;
		return erg;
	}

	////////////////////////////////////////////////////////////////////

	inline const value_type &entry(size_t i) const
	{
		return operator[] (i);
	}
	inline value_type &entry(size_t i)
	{
		return operator[] (i);
	}

	template<typename Type>
	DenseVector<TStorage> &
	assign(const Type &t);

	void maple_print(const char *name);
	
	
	inline size_t num_rows() const
	{
		return size();
	}
	
	inline size_t num_cols() const
	{
		return 1;
	}
	
	inline const value_type &operator() (size_t r, size_t c) const
	{
		UG_ASSERT(c==0, "vector only has one column");
		return entry(r);
	}
	
	inline value_type &operator() (size_t r, size_t c)
	{
		UG_ASSERT(c==0, "vector only has one column");
		return entry(r);
	}
	
	void
	subassign(size_t i, const value_type &t)
	{
		entry(i) = t;
	}
	
	template<typename T2>
	void
	subassign(size_t i, const T2 &t)
	{
		UG_ASSERT(i+t.size() <= size(), "");
		for(size_t i1=0; i1<t.size(); i1++)
			entry(i+i1) = t[i1];
	}
};

template<typename TStorage> 
inline
DenseVector<TStorage> 
operator * (double alpha, const DenseVector<TStorage> &vec)
{
	return vec*alpha;
}

template<typename TStorage>
inline
bool
operator > (double alpha, const DenseVector<TStorage> &vec)
{
	return vec < alpha;
}
template<typename TStorage>
inline
bool
operator < (double alpha, const DenseVector<TStorage> &vec)
{
	return vec > alpha;
}

template<typename TStorage>
inline
bool
operator >= (double alpha, const DenseVector<TStorage> &vec)
{
	return vec <= alpha;
}

template<typename TStorage>
inline
bool
operator <= (double alpha, const DenseVector<TStorage> &vec)
{
	return vec >= alpha;
}

// end group small_algebra
/// \}

}
#include "densevector_impl.h"

#endif // __H__UG__COMMON__DENSEVECTOR_H__

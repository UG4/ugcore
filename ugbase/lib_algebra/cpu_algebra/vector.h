/*
 *  Vector.h
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__CPU_ALGEBRA__VECTOR__
#define __H__UG__CPU_ALGEBRA__VECTOR__

#include "sparsematrix.h"

#include "../common/template_expressions.h"
#include "../common/operations.h"
#include "common/util/smart_pointer.h"
#include <vector>

namespace ug{
///////////////////////////////////////////////////////////////////
//							Vector
///////////////////////////////////////////////////////////////////

/// \addtogroup lib_algebra
///	@{

//!
template <typename TValueType>
class Vector
{
public:
	typedef TValueType value_type;
	//typedef subvector<value_type> subvector_type;
	typedef Vector<TValueType> vector_type;

	//! constructor
	Vector();

	//! constructor with size
	Vector(size_t size);

	//! virtual destructor
	virtual ~Vector();

	Vector(const vector_type & v)
	{
		m_size = 0; values = NULL;
		create(v.m_size);
		operator =(v);
	}

public:
	//! create a vector with specific size
	bool create(size_t size);
	//! create as a copy of other vector
	bool create(const Vector &v);
	
	//! clones the vector (deep-copy) including values
	SmartPtr<vector_type> clone() const;

	//! clones the vector (deep-copy) excluding values
	SmartPtr<vector_type> clone_without_values() const;

	//! resize vector.
	//! if bigger than capacity, new capacity is newSize+oldSize/2
	//! so growth factor is 1.5 (this is to keep memory overhead small)
	//! \see also https://github.com/facebook/folly/blob/master/folly/docs/FBVector.md
	void resize_sloppy(size_t newSize, bool bCopyValues=true);	
	
	//! resize the vector to be EXACTLY newSize big (no overhead)
	void resize_exactly(size_t newSize, bool bCopyValues=true);	
	
	//! reserve will allocate EXACTLY newCapacity
	//! assertion if newCapacity < size
	void reserve_exactly(size_t newCapacity, bool bCopyValues);
	
	//! reserve capacity in vector.
	//! if bigger than capacity, new capacity is newCapacity+oldCapacity/2
	void reserve_sloppy(size_t newCapacity, bool bCopyValues=true);
	
	void resize(size_t newSize, bool bCopyValues=true)
	{
		resize_exactly(newSize, bCopyValues);
	}
	void reserve(size_t newCapacity, bool bCopyValues=true)
	{
		reserve_exactly(newCapacity, bCopyValues);
	}
	
	
	size_t capacity()
	{
		return m_capacity;
	}

	//! access element i of the vector
	inline value_type &operator [] (size_t i);
	inline const value_type &operator [] (size_t i) const;


	//! returns v.T w, that is the dotprod of this vector and w
	double dotprod(const Vector &w); //const;

	// deprecated, use x.T() * y.
	//inline double operator *(const Vector &w); ///< shortcut for .dotprod(w)

	//double energynorm2(const SparseMatrix &A) const;
	/*double energynorm(const SparseMatrix &A) const
	{
		return sqrt(energynorm2(A));
	}*/

	//! assign double d to whole Vector
	double operator = (double d);
	//! assign double d to whole Vector
	bool set(double d) { operator = (d); return true; }
	bool set_random(double from, double to);

	/** add/set/get a local vector
	 *
	 * The local vector type must provide the following members:
	 * - size()					- length of local vector
	 * - index(size_t i)		- global index for component i
	 * - operator[](size_t i)	- access to value of component i
	 */
	template <typename V> bool add(const V& u);
	template <typename V> bool set(const V& u);
	template <typename V> bool get(V& u) const;


	bool add(const value_type *u, const size_t *indices, int nr);
	bool set(const value_type *u, const size_t *indices, int nr);
	bool get(value_type *u, const size_t *indices, int nr) const;


	//template<typename T> inline void apply(Operation_type op, const T &t);

	//! assign other vector v
	inline void operator = (const Vector &v);
	inline void operator += (const Vector &v);
	inline void operator -= (const Vector &v);

	inline bool operator *= (const number &a)
	{
		for(size_t i=0; i<size(); i++) values[i] *= a;
		return true;
	}

	//! return sqrt(sum values[i]^2) (euclidian norm)
	inline double norm() const;

	size_t size() const { return m_size; }


public: // output functions
	//! print vector to console
	void print(const char * const text = NULL) const;
	void p() {print(); } ///< gdb shortcut for print

	//! ostream << operator
	friend std::ostream &operator<<(std::ostream &output, const Vector &v)
	{
		output << "Vector " <<  "[" << v.m_size << "]";
		return output;
	}

	size_t defragment() { return true; }

public:
	/*size_t begin_index() { return 0;}
	size_t end_index() { return size();}

	value_type *begin() { return values + begin_index(); }
	value_type *end() { return values + end_index(); }*/

protected:
	//! virtual clone using covariant return type
	/**
	 * This should be used instead of the copy constructor at the places
	 * where the additional information stored in the object of the derived
	 * class (like the geometry or the topology of the grid) should be kept.
	 */
	virtual vector_type* virtual_clone() const;

	//! virtual clone using covariant return type excluding values
	/**
	 * This should be used instead of the copy constructor at the places
	 * where the additional information stored in the object of the derived
	 * class (like the geometry or the topology of the grid) should be kept.
	 */
	virtual vector_type* virtual_clone_without_values() const;

private:
	bool destroy();

	size_t m_size;			///< size of the vector (vector is from 0..size-1)
	size_t m_capacity;		///< size of the vector (vector is from 0..size-1)
	value_type *values;		///< array where the values are stored, size m_size

	//mutable vector_mode dist_mode;
};

// @}

} // namespace ug

#include "vector_impl.h"

#endif

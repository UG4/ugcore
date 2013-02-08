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

	//! constructor with length
	Vector(size_t _length);

	//! destructor
	virtual ~Vector();

	Vector(const vector_type & v)
	{
		length = 0; values = NULL;
		create(v.length);
		operator =(v);
	}

public:
	//! create a vector with specific length
	bool create(size_t _length);
	//! create as a copy of other vector
	bool create(const Vector &v);
	
	//! clone the vector (i.e. create a new one basing on the derived class)
	/**
	 * This should be used instead of the copy constructor at the places
	 * where the additional information stored in the object of the derived
	 * class (like the geometry or the topology of the grid) should be kept.
	 */
	//virtual SmartPtr<Vector> virtual_clone();

	//! resize vector
	bool resize(size_t new_length, bool bCopyValues=true);

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

	//! printofile: posx posy value
	void printtofile(const char *filename);

	size_t size() const { return length; }


public: // output functions
	//! print vector to console
	void print(const char * const text = NULL) const;
	void p() {print(); } ///< gdb shortcut for print

	//! ostream << operator
	friend std::ostream &operator<<(std::ostream &output, const Vector &v)
	{
		output << "Vector " <<  "[" << v.length << "]";
		return output;
	}

	size_t defragment() { return true; }

public:
	/*size_t begin_index() { return 0;}
	size_t end_index() { return size();}

	value_type *begin() { return values + begin_index(); }
	value_type *end() { return values + end_index(); }*/


private:
	bool destroy();

	size_t length;				///< length of the vector (vector is from 0..length-1)
	value_type *values;		///< array where the values are stored, size length

	//mutable vector_mode dist_mode;
};

// @}

} // namespace ug

#include "vector_impl.h"

#endif

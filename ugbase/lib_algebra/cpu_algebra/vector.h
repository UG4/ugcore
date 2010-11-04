/*
 *  Vector.h
 *  flexamg
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

namespace ug{
///////////////////////////////////////////////////////////////////
//							Vector
///////////////////////////////////////////////////////////////////

/// \addtogroup lib_algebra
///	@{

//!
//! "big" Vector class for use with the big SparseMatrix
//! can = template expressions like x = 0.5*x - y + A*z
//! see TemplateExpressions.h
template <typename templ_value_type>
class Vector :  //public TE_VEC<Vector<templ_value_type> >,
				public virtual IFunctionBase
{
public:
	typedef templ_value_type value_type;
	//typedef subvector<value_type> subvector_type;
	typedef Vector<templ_value_type> vector_type;

	//! constructor
	Vector();

	//! constructor with length
	Vector(size_t _length);

	//! destructor
	~Vector();

private: // forbidden functions
	Vector(Vector&); // disallow copy operator

public:
	//! create a vector with specific length
	bool create(size_t _length);
	//! create as a copy of other vector
	bool create(const Vector &v);
	//! resize vector
	bool resize(size_t new_size)
	{
		if(size() == new_size) return true;
		destroy();
		return create(new_size);
	}

	bool destroy();


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


	void add(const value_type &d, size_t i);
	void set(const value_type &d, size_t i);
	void get(value_type &d, size_t i) const;

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
	inline bool operator *= (const number &a)
	{
		for(size_t i=0; i<size(); i++) values[i] *= a;
		return true;
	}

	//! assign this vector to another vector v
	inline void applyto(Vector &v) const;

	// template expressions for functions
	template<class Function> inline void operator = (Function &ex);

	//! template expression assignment
	template<typename Type> inline void operator = (const Type &t);
	//! template expression +=
	template<typename Type> inline void operator += (const Type &t);
	//! template expression -=
	template<typename Type> inline void operator -= (const Type &t);

	//! return sqrt(sum values[i]^2) (euclidian norm)
	inline double norm();
	inline double two_norm() { return norm(); }

	//! printofile: posx posy value
	void printtofile(const char *filename);

	size_t size() const { return length; }

	void addTo(value_type &dest, size_t i) const
	{
		dest += values[i];
	}

	void substractFrom(value_type &dest, size_t i) const
	{
		dest -= values[i];
	}

	void assign(value_type &dest, size_t i) const
	{
		dest = values[i];
	}

	void preventForbiddenDestination(void *p, bool &bFirst) const
	{
		assert(bFirst == true || p != this);
		bFirst = false;
	}

public: // output functions
	//! print vector to console
	void print(const char * const text = NULL) const;
	void p() {print(); } ///< gdb shortcut for print

	//! ostream << operator
	friend ostream &operator<<(ostream &output, const Vector &v)
	{
		output << "Vector " <<  "[" << v.length << "]";
		return output;
	}

	void printtype() const;
	size_t finalize() { return true; }

public:
	int firstIndex() { return 0; }
	int lastIndex() { return size(); }


private:
	size_t length;				///< length of the vector (vector is from 0..length-1)
	value_type *values;		///< array where the values are stored, size length

	//mutable vector_mode dist_mode;
};

// @}

} // namespace ug

#include "vector_impl.h"

#endif

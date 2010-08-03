/*
 *  Vector.hpp
 *  flexamg
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__VECTOR_IMPL__
#define __H__UG__MARTIN_ALGEBRA__VECTOR_IMPL__

#include <fstream>
#include "algebra_misc.h"

#define prefetchReadWrite(a)

namespace ug{
template<typename entry_type>
inline entry_type &Vector<entry_type>::operator [] (size_t i)
{
	UG_ASSERT(i >= 0 && i < length, *this << ": tried to access element " << i);
	return values[i];
}

template<typename entry_type>
inline const entry_type &Vector<entry_type>::operator [] (size_t i) const
{
	UG_ASSERT(i >= 0 && i < length, *this << ": tried to access element " << i);
	return values[i];
}


// energynorm2 = x*(A*x)
/*inline double Vector<entry_type>::energynorm2(const SparseMatrix &A) const
{
	double sum=0;
	for(size_t i=0; i<length; i++)	sum += (A[i] * (*this)) * values[i];
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, sum += A[i] * (*this) * values[i]);
	return sum;
}*/

// dotprod
template<typename entry_type>
inline double Vector<entry_type>::dotprod(const Vector &w) //const
{
	UG_ASSERT(length == w.length,  *this << " has not same length as " << w);

	double sum=0;
	for(size_t i=0; i<length; i++)	sum += values[i] * w[i];
	return sum;
}

// assign double to whole Vector
template<typename entry_type>
inline double Vector<entry_type>::operator = (double d)
{
	for(size_t i=0; i<length; i++)
		values[i] = d;
	return d;
}



// für Function Expression, sh. TemplateExpression.h
template<typename entry_type>
template<class Function> inline void Vector<entry_type>::operator = (Function &ex)
{
	ex.applyto(*this);
}

template<typename entry_type>
inline void Vector<entry_type>::operator = (const Vector &v)
{
	v.applyto(*this);
}
template<typename entry_type>
inline void Vector<entry_type>::applyto(Vector &v) const
{
	UG_ASSERT(v.length == length, *this << " has not same length as " << v);

	for(size_t i=0; i<length; i++)
		v.values[i] = values[i];
}


template<typename entry_type>
template<typename Type> inline void Vector<entry_type>::operator = (const Type &t)
{
	//IF_PRINTLEVEL(5) cout << *this << " = " << t << " (unspecialized) " << endl;
	UG_ASSERT(t.size() == length, *this << " has not same length as " << t);
	t.preventForbiddenDestination(this);

	for(size_t i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		t.assign(values[i], i);
	}
}

// v += exp
template<typename entry_type>
template<typename Type> inline void Vector<entry_type>::operator += (const Type &t)
{
	UG_ASSERT(t.size() == length, *this << " has not same length as " << t);

	for(size_t i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		t.addTo(values[i], i);
	}
}

// v -= exp
template<typename entry_type>
template<typename Type> inline void Vector<entry_type>::operator -= (const Type &t)
{
	UG_DLOG(LIB_ALG_VECTOR, 5, *this << " -= " << t << " (unspecialized) ");
	UG_ASSERT(t.size() == length, *this << " has not same length as " << t);
	//t.preventForbiddenDestination(this);

	for(size_t i=0; i < length; i++)
		t.substractFrom(values[i], i);
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename entry_type>
Vector<entry_type>::Vector ()
{
	FORCE_CREATION { p(); } // force creation of this rountines for gdb.

	length = 0; values = NULL;
}

template<typename entry_type>
Vector<entry_type>::Vector(size_t _length)
{
	FORCE_CREATION { p(); } // force creation of this rountines for gdb.

	length = 0;
	create(_length);
}

template<typename entry_type>
Vector<entry_type>::~Vector()
{
	destroy();
}

template<typename entry_type>
bool Vector<entry_type>::destroy()
{
	if(values)
	{
		delete [] values;
		values = NULL;
	}
	length = 0;
	return true;
}


template<typename entry_type>
bool Vector<entry_type>::create(size_t _length)
{
	UG_ASSERT(length == 0, *this << " already created");
	length = _length;
	values = new entry_type[length];

	return true;
}


template<typename entry_type>
bool Vector<entry_type>::create(const Vector &v)
{
	UG_ASSERT(length == 0, *this << " already created");
	length = v.length;
	values = new entry_type[length];

	// we cannot use memcpy here bcs of variable blocks.
	for(size_t i=0; i<length; i++)
		values[i] = v.values[i];

	return true;
}


// print
template<typename entry_type>
void Vector<entry_type>::print(const char * const text) const
{

	if(text) cout << " == " << text;
	cout << " length: " << length << " =================" << endl;
	for(size_t i=0; i<length; i++)
		//cout << values[i] << " ";
		cout << i << ": " << values[i] << endl;
	cout << endl;
}

template<typename entry_type>
void Vector<entry_type>::printtype() const
{
	cout << *this;
}


template<typename entry_type>
template<typename V>
bool Vector<entry_type>::add(const V& u)
{
	for(size_t i=0; i < u.size(); i++)
		values[u.index(i)] += u[i];
	return true;
}

template<typename entry_type>
template<typename V>
bool Vector<entry_type>::set(const V& u)
{
	for(size_t i=0; i < u.size(); i++)
		values[u.index(i)] = u[i];
	return true;
}

template<typename entry_type>
template<typename V>
bool Vector<entry_type>::get(V& u) const
{
	for(size_t i=0; i < u.size(); i++)
		u[i] = values[u.index(i)];
	return true;
}




template<typename entry_type>
bool Vector<entry_type>::add(const entry_type *u, const size_t *indices, int nr)
{
	for(size_t i=0; i < nr; i++)
		values[indices[i]] += u[i];
	return true;
}

template<typename entry_type>
bool Vector<entry_type>::set(const entry_type *u, const size_t *indices, int nr)
{
	for(size_t i=0; i < nr; i++)
		values[indices[i]] = u[i];
	return true;
}

template<typename entry_type>
bool Vector<entry_type>::get(entry_type *u, const size_t *indices, int nr) const
{
	for(size_t i=0; i < nr; i++)
		u[i] = values[indices[i]] ;
	return true;
}



template<typename entry_type>
void Vector<entry_type>::add(const entry_type &d, size_t i)
{
	values[i] += d;
}
template<typename entry_type>
void Vector<entry_type>::set(const entry_type &d, size_t i)
{
	values[i] = d;
}
template<typename entry_type>
void Vector<entry_type>::get(entry_type &d, size_t i) const
{
	d = values[i];
}


template<typename entry_type>
double operator *(const TRANSPOSED<Vector<entry_type> > &x, const Vector<entry_type> &y)
{
	return x.T().dotprod(y);
}

template<typename entry_type>
inline double Vector<entry_type>::norm()
{
	double d=0;
	for(size_t i=0; i<size(); ++i)
		d+=BlockNorm2(values[i]);
	return sqrt(d);
}

}//namespace ug

#endif

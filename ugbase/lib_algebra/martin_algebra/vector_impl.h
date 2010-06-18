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
inline double Vector<entry_type>::dotprod(const Vector &w) const
{
	UG_ASSERT(length == w.length,  *this << " has not same length as " << w);
	
	double sum=0;
	for(size_t i=0; i<length; i++)	sum += values[i] * w[i];
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, sum += values[i] * w[i]);
	return sum;
}

// assign double to whole Vector
template<typename entry_type>
inline double Vector<entry_type>::operator = (double d)
{
	for(size_t i=0; i<length; i++)
	{
		prefetchReadWrite(values+i+512);
		values[i] = d;
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = d);
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
	//memcpy(v.values, values, length*sizeof(entry_type));
	for(size_t i=0; i<length; i++)
		v.values[i] = values[i];
}


/*void operator = (const Expression<SparseMatrix, Multiply_Operator, Vector> ex)
 {
 ASSERT2(ex.size() == length, *this << " has not same length as " << ex);
 const matrix &m = ex.l;
 const Vector &r = ex.r;
 //for(size_t i=0; i < length; i++) values[i] = m[i]*r;
 FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = m[i]*r);
 }	*/

template<typename entry_type> 
template<typename Type> inline void Vector<entry_type>::operator = (const Type &t)
{ 
	//IF_PRINTLEVEL(5) cout << *this << " = " << t << " (unspecialized) " << endl;
	UG_ASSERT(t.size() == length, *this << " has not same length as " << t);
	t.preventForbiddenDestination(this);

	for(size_t i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		//values[i] = t[i];
		t.assign(values[i], i);
		//t.copyTo(values[i], i);
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] = t[i]);
}

// v += exp
template<typename entry_type> 
template<typename Type> inline void Vector<entry_type>::operator += (const Type &t)
{ 
	//IF_PRINTLEVEL(5) cout << *this << " += " << t << " (unspecialized) " << endl;
	UG_ASSERT(t.size() == length, *this << " has not same length as " << t);
	//t.preventForbiddenDestination(this);
	
	for(size_t i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		//values[i] += t[i]; 
		t.addTo(values[i], i);
	}
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] += t[i]);
} 

// v -= exp
template<typename entry_type> 
template<typename Type> inline void Vector<entry_type>::operator -= (const Type &t)
{
	UG_DLOG(LIB_ALG_VECTOR, 5, *this << " -= " << t << " (unspecialized) ");
	UG_ASSERT(t.size() == length, *this << " has not same length as " << t);
	//t.preventForbiddenDestination(this);
	
	for(size_t i=0; i < length; i++)
	{
		prefetchReadWrite(values+i+512);
		//values[i] -= t[i]; 
		t.substractFrom(values[i], i);
	}			
	//FOR_UNROLL_FWD(i, 0, length, UNROLL, values[i] -= t[i]);
}
/*template<typename entry_type>
template<typename T> inline void Vector<entry_type>::apply(Operation_type op, const T &t)
{
	if(op == OPERATION_SET)
		operator = (t);
	else if(op == OPERATION_ADD)
		operator += (t);
	if(op == OPERATION_SUB)
		operator -= (t);
}*/


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename entry_type>
Vector<entry_type>::Vector (const char *_name)
{
	FORCE_CREATION { p(); } // force creation of this rountines for gdb.
		
	length = 0; values = NULL; name = _name; level = 0;
}	

template<typename entry_type>
Vector<entry_type>::Vector(size_t _length, const char *_name)
{
	FORCE_CREATION { p(); } // force creation of this rountines for gdb.

	length = 0;
	create(_length);
	name = _name;
	level = 0;
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


// printofile: posx posy value
template<typename entry_type>
void Vector<entry_type>::printtofile(const char *filename)
{
/*	fstream fil(filename, ios::out);
	
	for(size_t i=0; i < length; i++)
	{
		postype pos = GetPosForIndex(i);
		fil << pos.x << " " << pos.y << " " << values[i] << endl;
	}*/
	
}

// print
template<typename entry_type>
void Vector<entry_type>::print(const char * const text) const
{
  
  if(name) cout << endl << "================ " << name;
	if(text) cout << " == " << text;
	cout << " == level: " << level << " length: " << length << " =================" << endl;
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

/*template<typename entry_type>
void Vector<entry_type>::add(const subvector<entry_type> &subvec)
{
	for(size_t i=0; i < subvec.getNr(); i++)
		values[subvec.getIndex(i)] += subvec(i);
}

template<typename entry_type>
void Vector<entry_type>::set(const subvector<entry_type> &subvec)
{
	for(size_t i=0; i < subvec.getNr(); i++)
		values[subvec.getIndex(i)] = subvec(i);	
}

template<typename entry_type>
void Vector<entry_type>::get(subvector<entry_type> &subvec) const
{
	for(size_t i=0; i < subvec.getNr(); i++)
		subvec(i) = values[subvec.getIndex(i)];
}*/

template<typename entry_type>
bool Vector<entry_type>::add(const local_vector_type &u, const local_index_type &ind)
{
	for(std::size_t i=0; i < ind.size(); i++)
		values[ind[i][0]] += u[i];
	return true;
}

template<typename entry_type>
bool Vector<entry_type>::set(const local_vector_type &u, const local_index_type &ind)
{
	for(std::size_t i=0; i < ind.size(); i++)
		values[ind[i][0]] = u[i];
	return true;
}

template<typename entry_type>
bool Vector<entry_type>::get(local_vector_type &u, const local_index_type &ind) const
{
	for(std::size_t i=0; i < ind.size(); i++)
		u[i] = values[ind[i][0]];
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
	
}//namespace ug

#endif

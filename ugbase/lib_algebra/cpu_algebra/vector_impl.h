/*
 *  Vector.hpp
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */

#ifndef __H__UG__CPU_ALGEBRA__VECTOR_IMPL__
#define __H__UG__CPU_ALGEBRA__VECTOR_IMPL__

#include <fstream>
#include <algorithm>
#include "algebra_misc.h"
#include "common/math/ugmath.h"
#include "vector.h" // for urand

#define prefetchReadWrite(a)

namespace ug{
template<typename value_type>
inline value_type &Vector<value_type>::operator [] (size_t i)
{
	UG_ASSERT(i < m_size, *this << ": tried to access element " << i);
	return values[i];
}

template<typename value_type>
inline const value_type &Vector<value_type>::operator [] (size_t i) const
{
	UG_ASSERT(i < m_size, *this << ": tried to access element " << i);
	return values[i];
}


// energynorm2 = x*(A*x)
/*inline double Vector<value_type>::energynorm2(const SparseMatrix &A) const
{
	double sum=0;
	for(size_t i=0; i<m_size; i++)	sum += (A[i] * (*this)) * values[i];
	//FOR_UNROLL_FWD(i, 0, m_size, UNROLL, sum += A[i] * (*this) * values[i]);
	return sum;
}*/

// dotprod
template<typename value_type>
inline double Vector<value_type>::dotprod(const Vector &w) //const
{
	UG_ASSERT(m_size == w.m_size,  *this << " has not same size as " << w);

	double sum=0;
	for(size_t i=0; i<m_size; i++)	sum += VecProd(values[i], w[i]);
	return sum;
}

// assign double to whole Vector
template<typename value_type>
inline double Vector<value_type>::operator = (double d)
{
	for(size_t i=0; i<m_size; i++)
		values[i] = d;
	return d;
}

template<typename value_type>
inline bool Vector<value_type>::set_random(double from, double to)
{
	for(size_t i=0; i<size(); i++)
		for(size_t j=0; j<GetSize(values[i]); j++)
			BlockRef(values[i], j) = urand(from, to);
	return true;
}


template<typename value_type>
inline void Vector<value_type>::operator = (const vector_type &v)
{
	resize(v.size());
	for(size_t i=0; i<m_size; i++)
		values[i] = v[i];
}

template<typename value_type>
inline void Vector<value_type>::operator += (const vector_type &v)
{
	UG_ASSERT(v.size() == size(), "vector sizes must match! (" << v.size() << " != " << size() << ")");
	for(size_t i=0; i<m_size; i++)
		values[i] += v[i];
}

template<typename value_type>
inline void Vector<value_type>::operator -= (const vector_type &v)
{
	UG_ASSERT(v.size() == size(), "vector sizes must match! (" << v.size() << " != " << size() << ")");
	for(size_t i=0; i<m_size; i++)
		values[i] -= v[i];
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


template<typename value_type>
Vector<value_type>::Vector () : m_size(0), m_capacity(0), values(NULL)
{
	FORCE_CREATION { p(); } // force creation of this rountines for gdb.
}

template<typename value_type>
Vector<value_type>::Vector(size_t size) : m_size(0), m_capacity(0), values(NULL)
{
	FORCE_CREATION { p(); } // force creation of this rountines for gdb.
	create(size);
}

template<typename value_type>
Vector<value_type>::~Vector()
{
	destroy();
}

template<typename value_type>
bool Vector<value_type>::destroy()
{
	if(values)
	{
		delete [] values;
		values = NULL;
	}
	m_size = 0;
	return true;
}


template<typename value_type>
bool Vector<value_type>::create(size_t size)
{
	if(m_size == size) return true;
	destroy();
	m_size = size;
	values = new value_type[size];
	m_capacity = size;

	return true;
}


template<typename value_type>
Vector<value_type>* Vector<value_type>::virtual_clone() const
{
	return new Vector<value_type>(*this);
}

template<typename value_type>
SmartPtr<Vector<value_type> > Vector<value_type>::clone() const
{
	return SmartPtr<Vector<value_type> >(this->virtual_clone());
}

template<typename value_type>
Vector<value_type>* Vector<value_type>::virtual_clone_without_values() const
{
	return new Vector<value_type>(this->m_size);
}

template<typename value_type>
SmartPtr<Vector<value_type> > Vector<value_type>::clone_without_values() const
{
	return SmartPtr<Vector<value_type> >(this->virtual_clone_without_values());
}

template<typename value_type>
void Vector<value_type>::reserve_exactly(size_t newCapacity, bool bCopyValues)
{
	UG_ASSERT(newCapacity >= m_size, "use resize, then reserve_exactly");
	value_type *new_values = new value_type[newCapacity];	
	// we cannot use memcpy here bcs of variable blocks.
	if(values != NULL && bCopyValues)
	{
		for(size_t i=0; i<m_size; i++)
			std::swap(new_values[i], values[i]);
		for(size_t i=m_size; i<newCapacity; i++)
			new_values[i] = 0.0;
	}
	if(values) delete [] values;
	values = new_values;
	m_capacity = newCapacity;
}


template<typename value_type>
void Vector<value_type>::reserve_sloppy(size_t newCapacity, bool bCopyValues)
{
	if(newCapacity <= m_capacity) return;
	reserve_exactly(newCapacity + m_capacity/2, bCopyValues);
}

template<typename value_type>
void Vector<value_type>::resize_sloppy(size_t newSize, bool bCopyValues)
{
	if(newSize > m_capacity)
	{
		size_t newCapacity = m_size/2 + newSize;
		reserve_exactly(newCapacity, true);
	}
	m_size = newSize;
}

template<typename value_type>
void Vector<value_type>::resize_exactly(size_t newSize, bool bCopyValues)
{
	if(newSize > m_capacity)
		reserve_exactly(newSize, true);
	m_size = newSize;
}



template<typename value_type>
bool Vector<value_type>::create(const Vector &v)
{
	UG_ASSERT(m_size == 0, *this << " already created");
	m_size = v.m_size;
	values = new value_type[m_size];
	m_capacity = m_size;

	// we cannot use memcpy here bcs of variable blocks.
	for(size_t i=0; i<m_size; i++)
		values[i] = v.values[i];

	return true;
}


// print
template<typename value_type>
void Vector<value_type>::print(const char * const text) const
{

	if(text) std::cout << " == " << text;
	std::cout << " size: " << m_size << " =================" << std::endl;
	for(size_t i=0; i<m_size; i++)
		//cout << values[i] << " ";
		std::cout << i << ": " << values[i] << std::endl;
	std::cout << std::endl;
}


template<typename value_type>
template<typename V>
bool Vector<value_type>::add(const V& u)
{
	for(size_t i=0; i < u.size(); i++)
		values[u.index(i)] += u[i];
	return true;
}

template<typename value_type>
template<typename V>
bool Vector<value_type>::set(const V& u)
{
	for(size_t i=0; i < u.size(); i++)
		values[u.index(i)] = u[i];
	return true;
}

template<typename value_type>
template<typename V>
bool Vector<value_type>::get(V& u) const
{
	for(size_t i=0; i < u.size(); i++)
		u[i] = values[u.index(i)];
	return true;
}




template<typename value_type>
bool Vector<value_type>::add(const value_type *u, const size_t *indices, size_t nr)
{
	for(size_t i=0; i < nr; i++)
		values[indices[i]] += u[i];
	return true;
}

template<typename value_type>
bool Vector<value_type>::set(const value_type *u, const size_t *indices, size_t nr)
{
	for(size_t i=0; i < nr; i++)
		values[indices[i]] = u[i];
	return true;
}

template<typename value_type>
bool Vector<value_type>::get(value_type *u, const size_t *indices, size_t nr) const
{
	for(size_t i=0; i < nr; i++)
		u[i] = values[indices[i]] ;
	return true;
}


template<typename value_type>
double operator *(const TRANSPOSED<Vector<value_type> > &x, const Vector<value_type> &y)
{
	return x.T().dotprod(y);
}

template<typename value_type>
inline double Vector<value_type>::norm() const
{
	double d=0;
	for(size_t i=0; i<size(); ++i)
		d+=BlockNorm2(values[i]);
	return sqrt(d);
}

template<typename TValueType>
bool CloneVector(Vector<TValueType> &dest, const Vector<TValueType>& src)
{
	dest.resize(src.size());
	return true;
}

template<typename TValueType, class TOStream>
void Serialize(TOStream &buf, const Vector<TValueType> &v)
{
	Serialize(buf, v.size());
	for(size_t i=0; i < v.size(); i++)
		Serialize(buf, v[i]);
}

template<typename TValueType, class TIStream>
void Deserialize(TIStream &buf, Vector<TValueType> &v)
{
	size_t s = Deserialize<size_t>(buf);
	v.resize(s);
	for(size_t i=0; i < s; i++)
		Deserialize(buf, v[i]);
}

}//namespace ug

#endif

/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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
inline void Vector<value_type>::set_random(double from, double to)
{
	for(size_t i=0; i<size(); i++)
		for(size_t j=0; j<GetSize(values[i]); j++)
			BlockRef(values[i], j) = urand(from, to);
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
void Vector<value_type>::destroy()
{
	if(values)
	{
		delete [] values;
		values = NULL;
	}
	m_size = 0;
}


template<typename value_type>
void Vector<value_type>::create(size_t size)
{
	if(m_size == size) return;
	destroy();
	m_size = size;
	values = new value_type[size];
	m_capacity = size;
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
void Vector<value_type>::create(const Vector &v)
{
	UG_ASSERT(m_size == 0, *this << " already created");
	m_size = v.m_size;
	values = new value_type[m_size];
	m_capacity = m_size;

	// we cannot use memcpy here bcs of variable blocks.
	for(size_t i=0; i<m_size; i++)
		values[i] = v.values[i];
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
void Vector<value_type>::add(const V& u)
{
	for(size_t i=0; i < u.size(); i++)
		values[u.index(i)] += u[i];
}

template<typename value_type>
template<typename V>
void Vector<value_type>::set(const V& u)
{
	for(size_t i=0; i < u.size(); i++)
		values[u.index(i)] = u[i];
}

template<typename value_type>
template<typename V>
void Vector<value_type>::get(V& u) const
{
	for(size_t i=0; i < u.size(); i++)
		u[i] = values[u.index(i)];
}




template<typename value_type>
void Vector<value_type>::add(const value_type *u, const size_t *indices, size_t nr)
{
	for(size_t i=0; i < nr; i++)
		values[indices[i]] += u[i];
}

template<typename value_type>
void Vector<value_type>::set(const value_type *u, const size_t *indices, size_t nr)
{
	for(size_t i=0; i < nr; i++)
		values[indices[i]] = u[i];
}

template<typename value_type>
void Vector<value_type>::get(value_type *u, const size_t *indices, size_t nr) const
{
	for(size_t i=0; i < nr; i++)
		u[i] = values[indices[i]] ;
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

template<typename value_type>
inline double Vector<value_type>::maxnorm() const
{
	double d=0;
	for(size_t i=0; i<size(); ++i)
		d = std::max(d, BlockMaxNorm(values[i]));
	return d;
}

template<typename TValueType>
void CloneVector(Vector<TValueType> &dest, const Vector<TValueType>& src)
{
	dest.resize(src.size());
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

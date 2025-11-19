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


#ifndef __H__UG__COMMON__DENSEVECTOR_IMPL_H__
#define __H__UG__COMMON__DENSEVECTOR_IMPL_H__

#include "common/serialization.h"

namespace ug{

template<typename TStorage>
DenseVector<TStorage>::DenseVector(double alpha)
{
	operator=(alpha);
}

template<typename TStorage>
DenseVector<TStorage>::DenseVector(const DenseVector<TStorage> &rhs) : TStorage(rhs)
{
}

// operations with vectors
template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator = (const DenseVector<TStorage> &rhs)
{
	if(this == &rhs) return *this;
	if(size() != rhs.size()) resize(rhs.size());

	for(size_type i=0; i<size(); i++)
		entry(i) = rhs[i];
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator+=(const DenseVector &rhs)
{
	for(size_type i=0; i<size(); i++)
		entry(i) += rhs[i];
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator-=(const DenseVector &rhs)
{
	for(size_type i=0; i<size(); i++)
		entry(i) -= rhs[i];
	return *this;
}


// operations with scalars
template<typename TStorage>
template<typename T>
DenseVector<TStorage> &
DenseVector<TStorage>::operator=(const T &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) = alpha;
	return *this;
}

template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator+=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) += alpha;
	return *this;
}

template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator-=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) -= alpha;
	return *this;
}

template<typename TStorage>
template<typename T>
DenseVector<TStorage> &
DenseVector<TStorage>::operator*=(const T &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) *= alpha;
	return *this;
}


template<typename TStorage>
DenseVector<TStorage> &
DenseVector<TStorage>::operator/=(const typename DenseVector<TStorage>::value_type &alpha)
{
	for(size_t i=0; i<size(); i++)
		entry(i) /= alpha;

	return *this;
}

template<typename TStorage>
void DenseVector<TStorage>::maple_print(const char *name)
{
	UG_LOG(name << " = vector([");
	for(size_t i=0; i<size(); ++i)
	{
		if(i > 0) UG_LOG(", ");
		UG_LOG(entry(i));
	}
	UG_LOG("]);\n");
}

// views
// methods


template<typename TStorage>
std::ostream &operator << (std::ostream &out, const DenseVector<TStorage> &vec)
{
	out << "[";
	for(size_t i=0; i<vec.size(); i++)
		out << " " << vec[i];
	out << " ] ";
//	out << "(" << vec.size() << ")";
	return out;
}



template<typename TStorage>
template<typename Type>
DenseVector<TStorage> &
DenseVector<TStorage>::assign(const Type &t)
{
	VectorAssign(*this, t);
	return *this;
}



template<size_t n> inline void Serialize(std::ostream &buff, const DenseVector<FixedArray1<number, n> > &vec)
{
	buff.write((char*)&vec, sizeof(vec));
}

template<size_t n> inline void Deserialize(std::istream &buff, DenseVector<FixedArray1<number, n> > &vec)
{
	buff.read((char*)&vec, sizeof(vec));
}


template<typename T> inline void Serialize(std::ostream &buff, const DenseVector<VariableArray1<T> > &vec)
{
	size_t s = vec.size();
	buff.write((char*)&s, sizeof(s));
	for(size_t i=0; i<s; i++)
		Serialize(buff, vec[i]);
}

template<typename T> inline void Deserialize(std::istream &buff, DenseVector<VariableArray1<double> > &vec)
{
	size_t s;
	buff.read((char*)&s, sizeof(s));
	vec.resize(s);
	for(size_t i=0; i<s; i++)
		Deserialize(buff, vec[i]);
}


template<typename T >
inline bool IsFiniteAndNotTooBig(const DenseVector<T> &v)
{
	for(size_t r=0; r<v.size(); r++)
		if(IsFiniteAndNotTooBig(v[r]) == false) return false;

	return true;
}


//MAKE_TEMPLATE_OPERATORS_VECTOR2(typename TStorage, DenseVector<TStorage>);

}

#endif
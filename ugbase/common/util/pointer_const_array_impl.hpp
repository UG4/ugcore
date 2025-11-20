/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__pointer_array_impl__
#define __H__UG__pointer_array_impl__

#include <cstring>
#include <algorithm>

namespace ug
{

template <typename TPtr>
PointerConstArray<TPtr>::
PointerConstArray() :
	m_array(nullptr), m_data(nullptr), m_size(0), m_capacity(0)
{
}

template <typename TPtr>
PointerConstArray<TPtr>::
PointerConstArray(const PointerConstArray& pa)
{
	assign_pointer_const_array(pa);
}

template <typename TPtr>
PointerConstArray<TPtr>::
~PointerConstArray()
{
//	only free local memory...
	if(m_data)
		delete[] m_data;
}

template <typename TPtr>
PointerConstArray<TPtr>& PointerConstArray<TPtr>::
operator = (const PointerConstArray& pa)
{
//	free memory if necessary
	if(m_data)
		delete[] m_data;

//	assign the new class
	assign_pointer_const_array(pa);

//	return reference to this
	return *this;
}

template <typename TPtr>
void PointerConstArray<TPtr>::
assign_pointer_const_array(const PointerConstArray& pa)
{
	m_size = pa.m_size;
//	check whether pa points to an external array or not
	if(pa.m_array != pa.m_data){
	//	we'll point to the same external array
		m_array = pa.m_array;
		m_data = nullptr;
		m_capacity = 0;
	}
	else{
	//	we'll create our own array and copy the contents
		if(m_size){
			m_data = new TPtr[m_size];
			memcpy(m_data, pa.m_data, sizeof(TPtr) * m_size);
		}else
			m_data = nullptr;
		m_capacity = m_size;
		m_array = m_data;
	}
}

template <typename TPtr>
inline size_t PointerConstArray<TPtr>::
size() const
{
	return m_size;
}

template <typename TPtr>
inline bool PointerConstArray<TPtr>::
empty() const
{
	return m_size == 0;
}

template <typename TPtr>
inline TPtr const PointerConstArray<TPtr>::
operator [] (size_t i) const
{
	UG_ASSERT(i < m_size, "bad index!");
	return m_array[i];
}

template <typename TPtr>
void PointerConstArray<TPtr>::
set_external_array(TPtr const *array, size_t size, bool bCopy)
{
	if(bCopy){
	//don't call resize!
		if(size > m_capacity){
			delete[] m_data;
			m_data = new TPtr[size];
			m_capacity = size;
		}

		if(size > 0)
			memcpy(m_data, array, sizeof(TPtr) * size);

		m_array = m_data;
		m_size = size;
	}
	else{
		m_array = array;
		m_size = size;
	}
}

template <typename TPtr>
void PointerConstArray<TPtr>::
reserve(size_t capacity)
{
//	only copy old data, if we're actually using it
	reserve(capacity, m_array == m_data);
}

template <typename TPtr>
void PointerConstArray<TPtr>::
reserve(size_t capacity, bool copyOldData)
{
	if(capacity > m_capacity){
		PtrArray newData = new TPtr[capacity];

		if(m_data && copyOldData)
			memcpy(newData, m_data, sizeof(TPtr) * m_capacity);

		m_capacity = capacity;

	//	update pointers to new memory location and free old data.
	//	if m_array still points to an external array, we won't update the pointer.
		if(m_array == m_data)
			m_array = newData;

		delete[] m_data;
		m_data = newData;
	}
}

template <typename TPtr>
void PointerConstArray<TPtr>::
clear()
{
	m_size = 0;
}

template <typename TPtr>
void PointerConstArray<TPtr>::
push_back(TPtr p)
{
	using namespace std;
//	first make sure, that m_array is contained in m_data
	if(m_array != m_data){
	//	since we push something back, we double the capacity now
		reserve(max<size_t>(1, m_size * 2), false);
		if(m_size > 0)
			memcpy(m_data, m_array, sizeof(TPtr) * m_size);
		m_array = m_data;
	}

//	now check whether there's space left in m_data
	if(m_size >= m_capacity)
		reserve(max<size_t>(1, m_capacity * 2), true);

//	assign the new entry
	m_data[m_size] = p;
	++m_size;
}

}//	end of namespace

#endif

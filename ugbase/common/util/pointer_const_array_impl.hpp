#ifndef __H__UG__pointer_array_impl__
#define __H__UG__pointer_array_impl__

#include <cstring>
#include <algorithm>

namespace ug
{

template <class TPtr>
PointerConstArray<TPtr>::
PointerConstArray() :
	m_array(NULL), m_data(NULL), m_size(0), m_capacity(0)
{
}

template <class TPtr>
PointerConstArray<TPtr>::
PointerConstArray(const PointerConstArray& pa)
{
	assign_pointer_const_array(pa);
}

template <class TPtr>
PointerConstArray<TPtr>::
~PointerConstArray()
{
//	only free local memory...
	if(m_data)
		delete[] m_data;
}

template <class TPtr>
PointerConstArray<TPtr>& PointerConstArray<TPtr>::
operator=(const PointerConstArray& pa)
{
//	free memory if necessary
	if(m_data)
		delete[] m_data;

//	assign the new class
	assign_pointer_const_array(pa);

//	return reference to this
	return *this;
}

template <class TPtr>
void PointerConstArray<TPtr>::
assign_pointer_const_array(const PointerConstArray& pa)
{
	m_size = pa.m_size;
//	check whether pa points to an external array or not
	if(pa.m_array != pa.m_data){
	//	we'll point to the same external array
		m_array = pa.m_array;
		m_data = NULL;
		m_capacity = 0;
	}
	else{
	//	we'll create our own array and copy the contents
		if(m_size){
			m_data = new TPtr[m_size];
			memcpy(m_data, pa.m_data, sizeof(TPtr) * m_size);
		}else
			m_data = NULL;
		m_capacity = m_size;
		m_array = m_data;
	}
}

template <class TPtr>
inline size_t PointerConstArray<TPtr>::
size() const
{
	return m_size;
}

template <class TPtr>
inline bool PointerConstArray<TPtr>::
empty() const
{
	return m_size == 0;
}

template <class TPtr>
inline TPtr const PointerConstArray<TPtr>::
operator[](size_t i) const
{
	UG_ASSERT(i < m_size, "bad index!");
	return m_array[i];
}

template <class TPtr>
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

template <class TPtr>
void PointerConstArray<TPtr>::
reserve(size_t capacity)
{
//	only copy old data, if we're actually using it
	reserve(capacity, m_array == m_data);
}

template <class TPtr>
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

template <class TPtr>
void PointerConstArray<TPtr>::
clear()
{
	m_size = 0;
}

template <class TPtr>
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

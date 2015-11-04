// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.05.2011 (m,d,y)

#ifndef __H__UG__binary_buffer_impl__
#define __H__UG__binary_buffer_impl__

#include <cassert>
#include <cstring>
#include "vector_util.h"

namespace ug
{

inline size_t BinaryBuffer::capacity() const
{
	return m_data.size();
}

inline size_t BinaryBuffer::read_pos() const
{
	return m_readPos;
}

inline size_t BinaryBuffer::write_pos() const
{
	return m_writePos;
}

inline void BinaryBuffer::read(char* buf, size_t size)
{
//	make sure that we only read valid data
	assert(m_readPos + size <= m_data.size());
	assert(m_readPos + size <= m_writePos);

//	copy the data
	memcpy(buf, &m_data.front() + m_readPos, size);

//	adjust read-pos
	m_readPos += size;
}

inline void BinaryBuffer::write(const char* buf, size_t size)
{
//	make sure that our data buffer is big enough
	if(m_writePos + size > m_data.size()){
	//	if the size of the data to be written exceeds the current
	//	data-buffers size, then we increase the data buffer by
	//	the size of the now written data.
	//	If not, then we'll increase the data-buffer by simply doubling
	//	its memory, to minimize reallocation costs (this is just a
	//	heuristic)
		if(size > m_data.size())
			m_data.resize(m_data.size() + size);
		else
			m_data.resize(m_data.size() * 2);
	}

//	copy the data
	memcpy(&m_data.front() + m_writePos, buf, size);

//	adjust write pos
	m_writePos += size;
}

inline char* BinaryBuffer::buffer()
{
	return GetDataPtr(m_data);
}

inline bool BinaryBuffer::eof()
{
	return m_readPos >= m_writePos;
}

}//	end of namespace

#endif

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.05.2011 (m,d,y)

#include "binary_buffer.h"

namespace ug
{

BinaryBuffer::BinaryBuffer() :
	m_readPos(0), m_writePos(0)
{
}

BinaryBuffer::BinaryBuffer(size_t bufSize) :
	m_data(bufSize), m_readPos(0), m_writePos(0)
{
}

void BinaryBuffer::clear()
{
//	PLEASE NOTE: There are classes and algorithms that rely on the
//	fact, that clear does not release the memory in m_data.
//	Please do not change this behavior - or severe performance
//	drawbacks will be the result. (Note that this behavior can
//	also be found in the stl::vector implementations).
	m_readPos = m_writePos = 0;
}

void BinaryBuffer::reserve(size_t newSize)
{
	if(newSize > m_data.size())
		m_data.resize(newSize);
}

void BinaryBuffer::set_read_pos(size_t pos)
{
	m_readPos = pos;
}

void BinaryBuffer::set_write_pos(size_t pos)
{
	m_writePos = pos;
}

}//	end of namespace

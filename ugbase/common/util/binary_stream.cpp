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

#include "binary_stream.h"
#include "../log.h"

using namespace std;

namespace ug{

BinaryStreamBuffer::BinaryStreamBuffer()
{
	reserve(128);
}

void BinaryStreamBuffer::clear()
{
	resize(0);
}

size_t BinaryStreamBuffer::size() const
{
	if(egptr() < pptr())
		return pptr() - buffer();
	return egptr() - eback();
}

void BinaryStreamBuffer::reserve(size_t newSize)
{
	size_t getPos = 0;
	size_t putPos = 0;
	size_t readEndPos = 0;
	if(!m_dataBuf.empty()){
		getPos = gptr() - eback();
		readEndPos = egptr() - eback();
		putPos = pptr() - buffer();
	}

	if(newSize > m_dataBuf.size())
		m_dataBuf.resize(newSize);

	setg(buffer(), buffer() + getPos, buffer() + readEndPos);
	setp(buffer() + putPos, end());
}

void BinaryStreamBuffer::resize(size_t newSize)
{
	size_t getPos = 0;
	size_t putPos = 0;
	if(!m_dataBuf.empty()){
		getPos = gptr() - eback();
		putPos = pptr() - buffer();
	}

	if(newSize > m_dataBuf.size())
		m_dataBuf.resize(newSize);

	if(getPos > newSize)
		getPos = newSize;
	if(putPos > newSize)
		putPos = newSize;

	setg(buffer(), buffer() + getPos, buffer() + newSize);
	setp(buffer() + putPos, end());
}

void BinaryStreamBuffer::reset()
{
	size_t readEndPos = egptr() - eback();

	setg(buffer(), buffer(), buffer() + readEndPos);
	setp(buffer(), end());
}

void BinaryStreamBuffer::write_jump(size_t jumpSize)
{
	gbump(jumpSize);
}

void BinaryStreamBuffer::read_jump(size_t jumpSize)
{
	pbump(jumpSize);
}

size_t BinaryStreamBuffer::get_read_pos() const
{
	return gptr() - eback();
}

BinaryStreamBuffer::int_type BinaryStreamBuffer::overflow(int_type c)
{
	if(c != traits_type::eof())
	{
		size_t curSize = size();
		reserve(curSize * 2);
		m_dataBuf[curSize] = c;
		setp(pptr() + 1, epptr());
	}
	return c;
}

BinaryStreamBuffer::int_type BinaryStreamBuffer::underflow()
{
	if(egptr() < pptr()){
	//	new data has been written since last reserve or last read operation
		setg(eback(), gptr(), pptr());
		if(gptr() < egptr()){
			return m_dataBuf[gptr() - eback()];
		}
	}

	return traits_type::eof();
}

}// end of namespace

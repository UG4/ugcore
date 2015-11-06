/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

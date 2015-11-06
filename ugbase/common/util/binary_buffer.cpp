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

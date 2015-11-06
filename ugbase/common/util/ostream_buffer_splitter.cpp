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

#include "ostream_buffer_splitter.h"

using namespace std;

namespace ug{

OStreamBufferSplitter::
OStreamBufferSplitter() : m_buf1(NULL), m_buf2(NULL)
{
	setp(m_buf, m_buf + BUF_SIZE);
}

OStreamBufferSplitter::
OStreamBufferSplitter(std::streambuf* buf1, std::streambuf* buf2)
{
	set_buffers(buf1, buf2);
}

OStreamBufferSplitter::
~OStreamBufferSplitter()
{
	flush();
}

void OStreamBufferSplitter::flush()
{
	if(m_buf1 && m_buf2 && (pptr() != pbase())){
	//	write whats left in the buffer to the split-buffers
		m_buf1->sputn(pbase(), pptr() - pbase());
		m_buf2->sputn(pbase(), pptr() - pbase());
	}
	setp(m_buf, m_buf + BUF_SIZE);
}

void OStreamBufferSplitter::
set_buffers(std::streambuf* buf1, std::streambuf* buf2)
{
	flush();

	m_buf1 = buf1;
	m_buf2 = buf2;
}

ostream::int_type OStreamBufferSplitter::
overflow(ostream::int_type c)
{
	if(m_buf1 && m_buf2){
		m_buf1->sputn(pbase(), pptr() - pbase());
		m_buf2->sputn(pbase(), pptr() - pbase());
		if(!traits_type::eq_int_type(m_buf1->sputc(c), traits_type::eof())
			&& !traits_type::eq_int_type(m_buf2->sputc(c), traits_type::eof()))
		{
			setp(m_buf, m_buf + BUF_SIZE);
			return traits_type::not_eof(c);
		}
	}
	return traits_type::eof();
}

}// end of namespace

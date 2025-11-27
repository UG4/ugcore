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

#ifndef __H__UG__empty_stream__
#define __H__UG__empty_stream__

#include <iostream>

namespace ug
{

/// \addtogroup ugbase_common_io
/// \{

///	Used by EmptyStream, to send tokens into nirvana!
/**	Note that the implementation is not really the fastest, since a virtual
 * function is called for each token. However, not much is done in that
 * virtual function...
 */
class EmptyStreamBuffer : public std::streambuf
{
	public:
		EmptyStreamBuffer()
			{setp(m_buf, m_buf + BUF_SIZE);}

	//	implementation of virtual std::streambuf methods.
		inline int_type overflow(int_type c = traits_type::eof()) override {
			setp(m_buf, m_buf + BUF_SIZE);
			return c;
		}

	protected:
		static constexpr int BUF_SIZE = 128;
		char_type m_buf[BUF_SIZE];
};


///	a specialization of std::ostream, which doesn't write anything
/**	The intention of the empty-buffer is, that methods and functions in
 * a stream are executed, even if nothing shall be printed or be written to
 * a file. This is important since it allows us to to call communicating
 * methods during parallel output without risking a stall.
 */
class EmptyOStream : public std::ostream
{
	public:
		EmptyOStream() : std::ostream(&m_streamBuf)	{}

	protected:
		EmptyStreamBuffer	m_streamBuf;
};

// end group ugbase_common_io
/// \}

}//	end of namespace

#endif

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

#ifndef __H__UG_OSTREAM_BUFFER_SPLITTER__
#define __H__UG_OSTREAM_BUFFER_SPLITTER__

#include <iostream>
#include <fstream>
#include "empty_stream.h"

namespace ug
{

/// \addtogroup ugbase_common_io
/// \{

///	forwards data written to this stream-buffer to other stream buffers.
/**	Make sure that the buffers are initialized and ready for writing,
 *	before you intend to write anything to the buffer.
 */
class OStreamBufferSplitter : public std::streambuf
{
	public:
		OStreamBufferSplitter();

		OStreamBufferSplitter(std::streambuf* buf1, std::streambuf* buf2);

		~OStreamBufferSplitter() override;

	//	flushes the local buffer into the associated buffers
		void flush();

		void set_buffers(std::streambuf* buf1, std::streambuf* buf2);

		int_type overflow(int_type c = traits_type::eof()) override;

	private:
		static constexpr int BUF_SIZE = 128;
		std::streambuf*	m_buf1;
		std::streambuf*	m_buf2;
		char_type m_buf[BUF_SIZE];
};

// end group ugbase_common_io
/// \}

}// end of namespace

#endif

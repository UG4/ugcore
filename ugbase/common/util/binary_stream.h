/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UTIL__BINARY_STREAM__
#define __H__UTIL__BINARY_STREAM__

#include <iostream>
#include <vector>
#include "vector_util.h"

namespace ug
{

/// \addtogroup ugbase_common_io
/// \{

////////////////////////////////////////////////////////////////////////
///	A special version of a std::streambuf, writes data directly into a buffer that is accessible at any time.
/** This buffer stores the data directly in a growing array.
 * You may receive a pointer to that array using get_buffer().
 * The size of the buffer is returned by get_buffer_size().
 *
 * Note that the pointer returned from get_buffer may be invalidated
 * after new insertions into the buffer.
 * Using resize() you may reserve a buffer-section. Input- and output-positions
 * won't be affected by a resize.
 */
class BinaryStreamBuffer : public std::streambuf
{
	public:
		BinaryStreamBuffer();

	/// clears the data and resets the read and write positions.
		void clear();

	///	Similar to resize, however, doesn't alter the read-area.
		void reserve(size_t newSize);

	/// resizes the data-buffer but does not alter read and write positions.
	/**	A resize extends the readable area.*/
		void resize(size_t newSize);

	/// set read- and write-positions to the start of the buffer.
		void reset();

	/// returns a pointer to the front of the buffer or nullptr if the buffer is empty.
		inline char_type* buffer()				{return reinterpret_cast<char_type*>(GetDataPtr(m_dataBuf));}

		inline const char_type* buffer() const	{return reinterpret_cast<const char_type*>(GetDataPtr(m_dataBuf));}

	/// returns the size of the buffer in bytes.
		size_t size() const;

	/// advances the write-pointer by jumpSize bytes
		void write_jump(size_t jumpSize);

	/// advances the read-pointer by jumpSize bytes
		void read_jump(size_t jumpSize);

	/// returns the read-position
		size_t get_read_pos() const;
			
	//	implementation of virtual std::streambuf methods.
		virtual int_type overflow(int_type c = traits_type::eof());

		virtual int_type underflow();

	protected:
	///	returns pointer to the first entry behind the buffer.
		inline char_type* end()					{return buffer() + m_dataBuf.size();}
		inline const char_type* end() const		{return buffer() + m_dataBuf.size();}

	protected:
		std::vector<unsigned char>	m_dataBuf;
};

////////////////////////////////////////////////////////////////////////
///	a specialzation of std::iostream, that uses a \sa BinaryStreamBuffer as buffer.
class BinaryStream : public std::iostream
{
	public:
		BinaryStream() : std::iostream(&m_streamBuf)	{}
		BinaryStream(size_t newSize) : std::iostream(&m_streamBuf)	{resize(newSize);}

		inline void clear()	//< clears the data and resets the read and write positions.
			{m_streamBuf.clear();}

		inline void resize(size_t newSize) //< resizes the data-buffer but does not alter read and write positions.
			{m_streamBuf.resize(newSize);}

		inline void reset() //< set read- and write-positions to the start of the buffer.
			{m_streamBuf.reset();}

		inline void* buffer() //< returns a pointer to the front of the buffer.
			{return m_streamBuf.buffer();}

		inline size_t size() const//< returns the size of the buffer in bytes.
			{return m_streamBuf.size();}

		inline void write_jump(size_t jumpSize) //< advances the write-pointer by jumpSize bytes
			{m_streamBuf.write_jump(jumpSize);}

		inline void read_jump(size_t jumpSize) //< advances the read-pointer by jumpSize bytes
			{m_streamBuf.read_jump(jumpSize);}

		inline size_t read_pos() const
			{return m_streamBuf.get_read_pos();}

	///	returns true if there is more data left to read.
	/** \todo this method should be removed. The normal stream-methdos should be used
	 *		instead (eof,...). However - those do not seem to work properly in the moment.
	 */
		inline bool can_read_more()
			{return m_streamBuf.get_read_pos() < size();}
	protected:
		BinaryStreamBuffer	m_streamBuf;
};

// end group ugbase_common_io
/// \}

}//	end of namespace

#endif

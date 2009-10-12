//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m08 d11

#ifndef __H__UTIL__BINARY_STREAM__
#define __H__UTIL__BINARY_STREAM__

#include <iostream>
#include <vector>

namespace ug
{
//TODO:	The performance of the buffer could most likely be improved,
//		if it would use the streambufs internal buffer.

////////////////////////////////////////////////////////////////////////
///	A special version of a std::streambuf, writes data directly into a buffer that is accessible at any time.
/**
 * This buffer stores the data directly in a growing array.
 * You may receive a pointer to that array using get_buffer().
 * The size of the buffer is returned by get_buffer_size().
 * Note that the pointer returned from get_buffer may be invalidated
 * after new insertions into the buffer.
 * Using resize() you may reserve a buffer-section. Input- and output-positions
 * won't be affected by a resize.
 */
class BinaryStreamBuffer : public std::streambuf
{
	public:
		BinaryStreamBuffer() : m_readPos(0), m_writePos(0)
			{setbuf(NULL, 0); setg(0,0,0); setp(0,0);}

		inline void clear()	///< clears the data and resets the read and write positions.
			{resize(0); reset();}

		inline void resize(int newSize) ///< resizes the data-buffer but does not alter read and write positions.
			{m_dataBuf.resize(newSize);}

		inline void reset() ///< set read- and write-positions to the start of the buffer.
			{m_readPos = 0; m_writePos = 0;}

		inline char* buffer() ///< returns a pointer to the front of the buffer.
			{return &m_dataBuf.front();}

		inline int size() ///< returns the size of the buffer in bytes.
			{return m_dataBuf.size();}

		inline void write_jump(int jumpSize) ///< advances the write-pointer by jumpSize bytes
			{m_writePos += jumpSize;}

		inline void read_jump(int jumpSize) ///< advances the read-pointer by jumpSize bytes
			{m_readPos += jumpSize;}

	//	implementation of virtual std::streambuf methods.
		inline virtual int_type overflow(int_type c = traits_type::eof())
		{
			if(c != traits_type::eof())
			{
				if(m_writePos < m_dataBuf.size())
					m_dataBuf[m_writePos] = c;
				else
					m_dataBuf.push_back(c);
				++m_writePos;
			}
			return c;
		}

		inline virtual int_type uflow()
		{
			char tmp = underflow();
			m_readPos++;
			return tmp;
		}

		inline virtual int_type underflow()
		{
			if(m_readPos < m_dataBuf.size())
				return m_dataBuf[m_readPos];
			return traits_type::eof();
		}

	protected:
		std::vector<char>	m_dataBuf;
		int m_readPos;
		int m_writePos;
};

////////////////////////////////////////////////////////////////////////
///	a specialzation of std::iostream, that uses a \sa BinaryStreamBuffer as buffer.
class BinaryStream : public std::iostream
{
	public:
		BinaryStream() : std::iostream(&m_streamBuf)	{}
		BinaryStream(int newSize) : std::iostream(&m_streamBuf)	{resize(newSize);}

		inline void clear()	///< clears the data and resets the read and write positions.
			{m_streamBuf.clear();}

		inline void resize(int newSize) ///< resizes the data-buffer but does not alter read and write positions.
			{m_streamBuf.resize(newSize);}

		inline void reset() ///< set read- and write-positions to the start of the buffer.
			{m_streamBuf.reset();}

		inline char* buffer() ///< returns a pointer to the front of the buffer.
			{return m_streamBuf.buffer();}

		inline int size() ///< returns the size of the buffer in bytes.
			{return m_streamBuf.size();}

		inline void write_jump(int jumpSize) ///< advances the write-pointer by jumpSize bytes
			{m_streamBuf.write_jump(jumpSize);}

		inline void read_jump(int jumpSize) ///< advances the read-pointer by jumpSize bytes
			{m_streamBuf.read_jump(jumpSize);}

	protected:
		BinaryStreamBuffer	m_streamBuf;
};

}//	end of namespace

#endif

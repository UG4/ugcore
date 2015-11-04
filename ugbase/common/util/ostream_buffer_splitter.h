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

		~OStreamBufferSplitter();

	//	flushes the local buffer into the associated buffers
		void flush();

		void set_buffers(std::streambuf* buf1, std::streambuf* buf2);

		virtual int_type overflow(int_type c = traits_type::eof());

	private:
		static const int		BUF_SIZE = 128;
		std::streambuf*	m_buf1;
		std::streambuf*	m_buf2;
		char_type		m_buf[BUF_SIZE];
};

// end group ugbase_common_io
/// \}

}// end of namespace

#endif

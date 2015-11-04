// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 20.03.2012 (m,d,y)

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
		inline virtual int_type overflow(int_type c = traits_type::eof())
		{
			setp(m_buf, m_buf + BUF_SIZE);
			return c;
		}

	protected:
		static const int	BUF_SIZE = 128;
		char_type	m_buf[BUF_SIZE];
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

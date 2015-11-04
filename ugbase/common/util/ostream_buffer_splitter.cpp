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

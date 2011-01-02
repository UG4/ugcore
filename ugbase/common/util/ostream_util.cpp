// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y11 m01 d02

#include "ostream_util.h"

using namespace std;

namespace ug{

OStreamBufferSplitter::
OStreamBufferSplitter() : m_buf1(NULL), m_buf2(NULL)
{

}

OStreamBufferSplitter::
OStreamBufferSplitter(std::streambuf* buf1, std::streambuf* buf2)
{
	//set_buffers(buf1, buf2);
}

OStreamBufferSplitter::
~OStreamBufferSplitter()
{
}

void OStreamBufferSplitter::
set_buffers(std::streambuf* buf1, std::streambuf* buf2)
{
	m_buf1 = buf1;
	m_buf2 = buf2;
}

ostream::int_type OStreamBufferSplitter::
overflow(ostream::int_type c)
{
	if(m_buf1 && m_buf2
		&& !traits_type::eq_int_type(m_buf1->sputc(c), traits_type::eof())
		&& !traits_type::eq_int_type(m_buf2->sputc(c), traits_type::eof()))
	{
		return traits_type::not_eof(c);
	}
	return traits_type::eof();
}


ostream::int_type OStreamBufferEmpty::
overflow(ostream::int_type v)
{
	return v;
}


}// end of namespace

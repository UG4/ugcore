#include "binary_stream.h"
#include "../log.h"

using namespace std;

namespace ug{

BinaryStreamBuffer::BinaryStreamBuffer()
{
	reserve(128);
}

void BinaryStreamBuffer::clear()
{
	resize(0);
}

size_t BinaryStreamBuffer::size() const
{
	if(egptr() < pptr())
		return pptr() - buffer();
	return egptr() - eback();
}

void BinaryStreamBuffer::reserve(size_t newSize)
{
	size_t getPos = 0;
	size_t putPos = 0;
	size_t readEndPos = 0;
	if(!m_dataBuf.empty()){
		getPos = gptr() - eback();
		readEndPos = egptr() - eback();
		putPos = pptr() - buffer();
	}

	if(newSize > m_dataBuf.size())
		m_dataBuf.resize(newSize);

	setg(buffer(), buffer() + getPos, buffer() + readEndPos);
	setp(buffer() + putPos, end());
}

void BinaryStreamBuffer::resize(size_t newSize)
{
	size_t getPos = 0;
	size_t putPos = 0;
	if(!m_dataBuf.empty()){
		getPos = gptr() - eback();
		putPos = pptr() - buffer();
	}

	if(newSize > m_dataBuf.size())
		m_dataBuf.resize(newSize);

	if(getPos > newSize)
		getPos = newSize;
	if(putPos > newSize)
		putPos = newSize;

	setg(buffer(), buffer() + getPos, buffer() + newSize);
	setp(buffer() + putPos, end());
}

void BinaryStreamBuffer::reset()
{
	size_t readEndPos = egptr() - eback();

	setg(buffer(), buffer(), buffer() + readEndPos);
	setp(buffer(), end());
}

void BinaryStreamBuffer::write_jump(size_t jumpSize)
{
	gbump(jumpSize);
}

void BinaryStreamBuffer::read_jump(size_t jumpSize)
{
	pbump(jumpSize);
}

size_t BinaryStreamBuffer::get_read_pos() const
{
	return gptr() - eback();
}

BinaryStreamBuffer::int_type BinaryStreamBuffer::overflow(int_type c)
{
	if(c != traits_type::eof())
	{
		size_t curSize = size();
		reserve(curSize * 2);
		m_dataBuf[curSize] = c;
		setp(pptr() + 1, epptr());
	}
	return c;
}

BinaryStreamBuffer::int_type BinaryStreamBuffer::underflow()
{
	if(egptr() < pptr()){
	//	new data has been written since last reserve or last read operation
		setg(eback(), gptr(), pptr());
		if(gptr() < egptr()){
			return m_dataBuf[gptr() - eback()];
		}
	}

	return traits_type::eof();
}

}// end of namespace

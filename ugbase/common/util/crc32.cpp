// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 18.02.2012 (d.m.y())

#include <cstring>
#include "boost/crc.hpp"
#include "crc32.h"

namespace ug{

uint32 crc32(const char* str)
{
	static boost::crc_32_type crc32Calculator;
	size_t len = strlen(str);
	crc32Calculator.process_bytes(str, len);
	uint32 checksum = crc32Calculator.checksum();
	crc32Calculator.reset();
	return checksum;
}

}//	end of namespace


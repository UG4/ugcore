// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 18.02.2012 (d.m.y)

#ifndef __H__UG_crc32__
#define __H__UG_crc32__

#include "../types.h"

namespace ug{

///	Calculates the crc32 for a null-terminated string.
uint32 crc32(const char* str);

}//	end of namespace

#endif


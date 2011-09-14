// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)

#ifndef __H__UG__os_info__
#define __H__UG__os_info__

namespace ug
{

///	returns the standard prefix of static and dynamic libraries on this os
const char* GetDynamicLibraryPrefix();

///	returns the standard suffix of dynamic libraries on this os
const char* GetDynamicLibrarySuffix();

}//	end of namespace

#endif

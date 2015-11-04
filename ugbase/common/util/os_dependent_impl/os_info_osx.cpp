// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)
 
#include "../os_info.h"

namespace ug{

const char* GetDynamicLibraryPrefix()
{
	return "lib";
}

const char* GetDynamicLibrarySuffix()
{
	return "dylib";
}

const char* GetPathSeparator()
{
	return "/";
}

}// end of namespace

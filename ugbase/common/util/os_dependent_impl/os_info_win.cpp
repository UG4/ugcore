// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)
 
#include "../os_info.h"

namespace ug{

const char* GetDynamicLibraryPrefix()
{
	return "";
}

const char* GetDynamicLibrarySuffix()
{
	return "dll";
}

const char* GetPathSeparator()
{
	return "\\";
}

}// end of namespace

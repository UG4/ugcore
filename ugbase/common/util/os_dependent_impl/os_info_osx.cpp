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

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

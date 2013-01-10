// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 10.01.2013 (d,m,y)

#include "../dynamic_library_util.h"

namespace ug{

DynLibHandle OpenLibrary(const char* fileName)
{
	return LoadLibrary(fileName);
}

bool CloseLibrary(DynLibHandle h)
{
	return FreeLibrary(h);
}

void* GetLibraryProcedure(DynLibHandle h, const char* procName)
{
	return GetProcAddress(h, procName);
}

}// end of namespace

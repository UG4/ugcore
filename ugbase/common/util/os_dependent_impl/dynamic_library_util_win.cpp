#include <string>
#include "../dynamic_library_util.h"

namespace ug{


/*std::string GetLastErrorString()
{
	DWORD errCode = GetLastError();
	char *err;
	if (!FormatMessage(FORMAT_MESSAGE_ALLOCATE_BUFFER | FORMAT_MESSAGE_FROM_SYSTEM | FORMAT_MESSAGE_IGNORE_INSERTS,
					   NULL,
					   errCode,
					   MAKELANGID(LANG_NEUTRAL, SUBLANG_DEFAULT), // default language
					   (LPTSTR) &err,
					   0,
					   NULL))
		return "?";

	std::string s = err;
	LocalFree(err);
	return s;
}
*/

DynLibHandle OpenLibrary(const char* fileName)
{
	DynLibHandle dlh = LoadLibrary(fileName);
	if(dlh == NULL)
		throw std::string("LoadLibrary failed."); //GetLastErrorString();
	return dlh;
}

bool CloseLibrary(DynLibHandle h)
{
	return FreeLibrary(h);
}

void* GetLibraryProcedure(DynLibHandle h, const char* procName)
{
	return (void*) GetProcAddress(h, procName);
}

}// end of namespace

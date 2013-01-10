// created by Sebastian Reiter
// s.b.reiter@gmail.com
// 10.01.2013 (d,m,y)

#ifndef __H__UG__dynamic_library_util__
#define __H__UG__dynamic_library_util__

namespace ug{

#ifdef UG_WIN32
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
///	Defines a reference to a dynamically loaded library
	typedef HMODULE DynLibHandle;
#else
///	Defines a reference to a dynamically loaded library
	typedef void*	DynLibHandle;
#endif

///	Loads a library and returns a handle
/**	The function returns 0 if no library with the specified name was found.
 * \sa CloseLibrary, GetLibraryProcedure*/
DynLibHandle OpenLibrary(const char* fileName);

///	Frees a library
/**	The method returns true if the operation was successful, false if not.
 * \sa OpenLibrary, GetLibraryProcedure*/
bool CloseLibrary(DynLibHandle h);

///	Returns the address of the specified procedure in the given library
/**	If no procedure with the given name was found, NULL is returned.
 * \sa OpenLibrary, CloseLibrary*/
void* GetLibraryProcedure(DynLibHandle h, const char* procName);

}// end of namespace

#endif

/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__dynamic_library_util__
#define __H__UG__dynamic_library_util__

#ifdef UG_WIN32
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
#endif

namespace ug{

/// \addtogroup ugbase_common_util
/// \{

#ifdef UG_WIN32
///	Defines a reference to a dynamically loaded library
	typedef HMODULE DynLibHandle;
#else
///	Defines a reference to a dynamically loaded library
	typedef void*	DynLibHandle;
#endif

///	Loads a library and returns a handle
/**	The function throws an std::string as error if the library could not be found
 * or if the library could not be opened (e.g. dependent shared libraries not found)
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

// end group ugbase_common_util
/// \}

}// end of namespace

#endif

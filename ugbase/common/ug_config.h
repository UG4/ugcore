/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__symbol_import_export__
#define __H__UG__symbol_import_export__

//	Those macros can be used if a class or function has to be
//	explicitly marked as a exported or imported method.
//	Currently this is only important on the windows platform, where they can
//	be used to appoint the library, which defines the function and the libraries
//	or executables, which just use the function.

/// \addtogroup ugbase_common
/// \{

#if defined UG_WIN32
	#ifdef __GNUC__
		#define EXPORT_IMPL __attribute__ ((dllexport))
		#define IMPORT_IMPL __attribute__ ((dllimport))
	#else
		#define EXPORT_IMPL __declspec(dllexport)
		#define IMPORT_IMPL __declspec(dllimport)
	#endif
#else
	#define EXPORT_IMPL
	#define IMPORT_IMPL
#endif


#ifdef BUILDING_DYNAMIC_LIBRARY
	#define UG_API EXPORT_IMPL
#else
	#ifdef IMPORTING_DYNAMIC_LIBRARY
		#define UG_API IMPORT_IMPL
	#else
		#define UG_API
	#endif
#endif

// end group ugbase_common
/// \}

#endif

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 19.09.2011 (m,d,y)

#ifndef __H__UG__symbol_import_export__
#define __H__UG__symbol_import_export__

//	Those macros can be used if a class or function has to be
//	explicitly marked as a exported or imported method.
//	Currently this is only important on the windows platform, where they can
//	be used to appoint the library, which defines the function and the libraries
//	or executables, which just use the function.


#if defined _WIN32 || defined __CYGWIN__
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
	#define UG_API IMPORT_IMPL
#endif

#endif

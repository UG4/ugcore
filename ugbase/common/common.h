//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#include <iostream>
#include <fstream>

#ifndef __H__LIB_GRID__COMMON__
#define __H__LIB_GRID__COMMON__

//	if LG_DEF__ENABLE_LOGGING is defined, all LOG(...) calls are written to libGrid_LOG.txt
//	if it is not defined, LOG calls are completely ignored (not even compiled).
#define LG_DEF__ENABLE_LOGGING

//	if we are compiling a library we want to log to a file.
//	define LG_DEF__COMPILE_LIBRARY in your makefile or your IDE
#ifdef LG_DEF__COMPILE_LIBRARY
	#define LG_DEF__LOG_TO_FILE
#endif

//	if LG_DEF__COMPILE_DEBUG is defined, special code will be compiled that performs checks on the validity of some operations.
//	Warning: if LG_DEF__COMPILE_DEBUG is defined, execution speed may suffer significantly.
#define LG_DEF__COMPILE_DEBUG

////////////////////////////////////////////////////////////////////////////////////////////////
//	LOG
///	if LG_DEF__ENABLE_LOGGING is defined, messages passed to LOG are directly written to libGrid_LOG.txt.
/**
 * The message passed to LOG can be of any type, that supports the iostream operator <<.
 * Multiple messages can be combined by connecting them with the << operator.
 * The following two lines show a valid application of LOG:
 * float year = 2008;
 * LOG("Development of libGrid started in " << year << " at the GCSC Frankfurt." << std::endl);
 */
#ifdef LG_DEF__ENABLE_LOGGING
	#ifdef LG_DEF__LOG_TO_FILE
		#include <fstream>
		std::ofstream& lib_grid_logger();
		#define LOG(msg) lib_grid_logger() << msg; lib_grid_logger().flush();
	#else
		#include <iostream>
		#define LOG(msg) std::cout << msg; std::cout.flush();
	#endif
#else
	#define LOG(msg)
#endif


////////////////////////////////////////////////////////////////////////////////////////////////
#define SAFE_DELETE(a)		{if(a){ delete a; a = NULL;}}
#define SAFE_RELEASE(p)		{if(p) { (p)->Release(); (p)=NULL;}}




#endif

//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d10

#ifndef __H__COMMON__COMMON__
#define __H__COMMON__COMMON__

#include <iostream>
#include <fstream>

/////////////////////////////////////////////////////////////////
// defines (currently here, should be compiling options)

// Warning and debug logs are disabled when compiling for release
#ifndef NDEBUG
	//#define UG_LOG_TO_FILE
	#define UG_ENABLE_WARNINGS
	#define UG_ENABLE_DEBUG_LOGS
#else /* NDEBUG */
	#undef UG_LOG_TO_FILE
	#undef UG_ENABLE_WARNINGS
	#undef UG_ENABLE_DEBUG_LOGS
#endif /* NDEBUG*/

/////////////////////////////////////////////////////////////////
// includes

#include "types.h"
//#include "contract.h"
#include "log.h"
#include "assert.h"
#include "static_assert.h"
#include "metaprogramming_util.h"

// depreciated, currently here for backward compatibility
#define LOG(msg) UG_LOG(msg)
#define STATIC_ASSERT(expr, msg) UG_STATIC_ASSERT(expr, msg)

////////////////////////////////////////////////////////////////////////////////////////////////
// save pointer (de-)allocation

#define SAFE_DELETE(a)		{if(a){ delete a; a = NULL;}}
#define SAFE_RELEASE(p)		{if(p) { (p)->Release(); (p)=NULL;}}

#endif

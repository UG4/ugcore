/*
The zlib/libpng License

Copyright (c) 2007 Aidin Abedi (www.*)

This software is provided 'as-is', without any express or implied warranty. In no event will
the authors be held liable for any damages arising from the use of this software.

Permission is granted to anyone to use this software for any purpose, including commercial 
applications, and to alter it and redistribute it freely, subject to the following
restrictions:

    1. The origin of this software must not be misrepresented; you must not claim that 
       you wrote the original software. If you use this software in a product, 
       an acknowledgment in the product documentation would be appreciated but is 
       not required.

    2. Altered source versions must be plainly marked as such, and must not be 
       misrepresented as being the original software.

    3. This notice may not be removed or altered from any source distribution.
*/

/*
This file has been altered by Sebastian Reiter (s.b.reiter@googlemail.com).
I replaced some includes and defines, to make it easier to integrate shiny
into ug.
New code is marked with //sreiter
*/

#ifndef SHINY_PREREQS_H
#define SHINY_PREREQS_H

#include "ShinyConfig.h"

#include "common/types.h"//sreiter
#include "stdlib.h"//sreiter

/*//sreiter
#if SHINY_PLATFORM == SHINY_PLATFORM_POSIX
#include <sys/types.h>
#endif
*///sreiter

namespace Shiny {


//-----------------------------------------------------------------------------
	
#if SHINY_PROFILER == TRUE
	struct ProfileNode;
	struct ProfileZone;

	using ProfileNodeCache = ProfileNode*;
	using ProfileNodeTable = ProfileNode*;
#endif


//-----------------------------------------------------------------------------


#if SHINY_COMPILER == SHINY_COMPILER_MSVC
#	define SHINY_INLINE		__forceinline
#	define SHINY_UNUSED		

#elif SHINY_PLATFORM == SHINY_COMPILER_GNUC
#	define SHINY_INLINE		__inline
#	define SHINY_UNUSED		__attribute__ ((unused))

//#elif SHINY_PLATFORM == SHINY_COMPILER_OTHER	//sreiter
#else
#	define SHINY_INLINE		inline
#	define SHINY_UNUSED		
#endif

//-----------------------------------------------------------------------------
/*//sreiter
#if SHINY_COMPILER == SHINY_COMPILER_MSVC
	using  int32_t = int;
	using uint32_t = unsigned int;

	using int64_t = __int64;
	using uint64_t = unsigned __int64;

#elif defined(__CYGWIN__)
	using uint32_t = u_int32_t;
	using uint64_t = u_int64_t;
#endif
*///sreiter
	using int32_t = int32;
	using uint32_t = uint32;
	using int64_t = int64;
	using uint64_t = uint64;

	using tick_t = uint64_t;

} // namespace Shiny

#endif // ifndef SHINY_*_H

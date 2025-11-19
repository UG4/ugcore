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
The hash-function problem is not yet solved!
The hash-function won't work correctly for 64-bit.
New code is marked with //sreiter
*/

#ifndef SHINY_TOOLS_H
#define SHINY_TOOLS_H

#include "ShinyPrereqs.h"


namespace Shiny {


//-----------------------------------------------------------------------------

	struct TimeUnit {
		float tickFreq;
		float invTickFreq;
		const char* suffix;
	};

	const TimeUnit* GetTimeUnit(float ticks);


//-----------------------------------------------------------------------------

	void GetTicks(tick_t *p);

	tick_t GetTickFreq(void);

	float GetTickInvFreq(void);


//-----------------------------------------------------------------------------

#if SHINY_COMPILER == SHINY_COMPILER_MSVC
#	pragma warning(disable: 4311)
#	pragma warning(disable: 4302)
#endif

	inline uint32_t ptr32(const void *a_Ptr) {
		unsigned long int tmp = reinterpret_cast<unsigned long int>(a_Ptr);//sreiter
		//uint32_t u = (uint32_t)tmp; //variable u has not been used, commented out by avogel, Dec 18 2009//sreiter
		return tmp;//sreiter
	}

#if SHINY_COMPILER == SHINY_COMPILER_MSVC
#	pragma warning(default: 4311)
#endif


} // namespace Shiny

#endif

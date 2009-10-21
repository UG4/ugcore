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

#include "ShinyTools.h"

#if SHINY_PLATFORM == SHINY_PLATFORM_WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>

#elif SHINY_PLATFORM == SHINY_PLATFORM_POSIX
#include <sys/time.h>
#endif

namespace Shiny {


//-----------------------------------------------------------------------------

	const TimeUnit* GetTimeUnit(float ticks) {
		static TimeUnit unit[] = {
			{ GetTickFreq() / 1.0f			, GetTickInvFreq() * 1.0f			, "s" },
			{ GetTickFreq() / 1000.0f		, GetTickInvFreq() * 1000.0f		, "ms" },
			{ GetTickFreq() / 1000000.0f	, GetTickInvFreq() * 1000000.0f		, "us" },
			{ GetTickFreq() / 1000000000.0f	, GetTickInvFreq() * 1000000000.0f	, "ns" }
		};

		if (unit[0].tickFreq < ticks) return &unit[0];
		else if (unit[1].tickFreq < ticks) return &unit[1];
		else if (unit[2].tickFreq < ticks) return &unit[2];
		else return &unit[3];
	}


//-----------------------------------------------------------------------------

#if SHINY_PLATFORM == SHINY_PLATFORM_WIN32

	tick_t _InitTickFreq(void) {
		tick_t freq;

		QueryPerformanceFrequency(reinterpret_cast<LARGE_INTEGER*>(&freq));
		return freq;
	}

	void GetTicks(tick_t *p) {
		QueryPerformanceCounter(reinterpret_cast<LARGE_INTEGER*>(p));
	}

	tick_t GetTickFreq(void) {
		static tick_t freq = _InitTickFreq();
		return freq;
	}

	float GetTickInvFreq(void) {
		static float invfreq = 1.0f / GetTickFreq();
		return invfreq;
	}


//-----------------------------------------------------------------------------

#elif SHINY_PLATFORM == SHINY_PLATFORM_POSIX

	void GetTicks(tick_t *p) {
		timeval time;
		gettimeofday(&time, NULL);

		*p = time.tv_sec * 1000000 + time.tv_usec;
	}

	tick_t GetTickFreq(void) {
		return 1000000;
	}

	float GetTickInvFreq(void) {
		return 1.0f / 1000000.0f;
	}

#endif
} // namespace Shiny

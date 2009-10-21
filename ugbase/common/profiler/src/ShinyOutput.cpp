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

#include "ShinyOutput.h"

#include <stdio.h>

#if SHINY_COMPILER == SHINY_COMPILER_MSVC
#	pragma warning(disable: 4996)
#	define snprintf		_snprintf
#	define TRAILING		0

#else
#	define TRAILING		1
#endif

#if SHINY_PROFILER == TRUE
namespace Shiny {


//-----------------------------------------------------------------------------

	void _printHeader(char *dest, const char *a_title) {
		snprintf(dest, OUTPUT_WIDTH_SUM + TRAILING,
			"%-*s %*s %*s %*s",
			OUTPUT_WIDTH_NAME, a_title,
			OUTPUT_WIDTH_HIT, "hits",
			OUTPUT_WIDTH_TIME+4+OUTPUT_WIDTH_PERC+1, "self time",
			OUTPUT_WIDTH_TIME+4+OUTPUT_WIDTH_PERC+1, "total time");
	}


//-----------------------------------------------------------------------------

	void _printData(char *dest, const ProfileData &a_data, float a_topercent) {
		float totalTicksAvg = a_data.totalTicksAvg();
		const TimeUnit *selfUnit = GetTimeUnit(a_data.selfTicks.avg);
		const TimeUnit *totalUnit = GetTimeUnit(totalTicksAvg);

		snprintf(dest, OUTPUT_WIDTH_DATA + TRAILING,
			" %*.1f %*.0f %-2s %*.0f%% %*.0f %-2s %*.0f%%",
			OUTPUT_WIDTH_HIT, a_data.entryCount.avg,
			OUTPUT_WIDTH_TIME, a_data.selfTicks.avg * selfUnit->invTickFreq, selfUnit->suffix,
			OUTPUT_WIDTH_PERC, a_data.selfTicks.avg * a_topercent,
			OUTPUT_WIDTH_TIME, totalTicksAvg * totalUnit->invTickFreq, totalUnit->suffix,
			OUTPUT_WIDTH_PERC, totalTicksAvg * a_topercent);
	}

//-----------------------------------------------------------------------------

	std::string OutputNodesAsString(const ProfileNode *a_root, uint32_t a_count) {
		float fTicksToPc = 100.0f / a_root->data.childTicks.avg;
		std::string str;

		str.resize((1 + a_count) * (OUTPUT_WIDTH_SUM + 1) - 1);
		char *s = &str[0];

		_printHeader(s, "call tree");
		s += OUTPUT_WIDTH_SUM;
		(*s++) = '\n';

		const ProfileNode *node = a_root;

		do {
			int offset = node->entryLevel * 2;
			snprintf(s, OUTPUT_WIDTH_NAME + TRAILING, "%*s%-*s",
				offset, "", OUTPUT_WIDTH_NAME - offset, node->zone->name);

			s += OUTPUT_WIDTH_NAME;

			_printData(s, node->data, fTicksToPc);

			s += OUTPUT_WIDTH_DATA;
			(*s++) = '\n';

			node = node->findNextInTree();
		} while (node);

		*(--s) = '\0';
		return str;
	}


//-----------------------------------------------------------------------------

	std::string OutputZonesAsString(const ProfileZone *a_root, uint32_t a_count) {
		float fTicksToPc = 100.0f / a_root->data.childTicks.avg;
		std::string str;

		str.resize((1 + a_count) * (OUTPUT_WIDTH_SUM + 1) - 1);
		char *s = &str[0];

		_printHeader(s, "flat profile");
		s += OUTPUT_WIDTH_SUM;
		(*s++) = '\n';

		const ProfileZone *zone = a_root;

		do {
			snprintf(s, OUTPUT_WIDTH_NAME + TRAILING, "%-*s",
				OUTPUT_WIDTH_NAME, zone->name);

			s += OUTPUT_WIDTH_NAME;

			_printData(s, zone->data, fTicksToPc);

			s += OUTPUT_WIDTH_DATA;
			(*s++) = '\n';

			zone = zone->next;
		} while (zone);

		*(--s) = '\0';
		return str;
	}

} // namespace Shiny
#endif

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

/* CHANGES!
 * This file contains changes concerning ProfileData::computeAverages. The old
 * damping (avg = a_damping * (avg - cur) + cur) was replaced by a new damping:
 * (avg = a_damping * avg + cur). The new daming allows to perform intermediate
 * updates while still conserving all gathered profile-times (choose a_damping==1). 
 * To only consider new profile-times (thus clearing the profile-history on update)
 * you may choose a_damping==0.
 * Please contact sreiter@gcsc.uni-frankfurt.de or mrupp@gcsc.uni-frankfurt.de
 * for more information and discussion on this subject.*/


#ifndef SHINY_DATA_H
#define SHINY_DATA_H

#include "ShinyPrereqs.h"

namespace Shiny {


//-----------------------------------------------------------------------------
	
	struct ProfileLastData {
		uint32_t entryCount;
		tick_t selfTicks;
	};


//-----------------------------------------------------------------------------

	struct ProfileData {

		template <typename T>
		struct Data {
			T cur;
			float avg;

		// CHANGE:	changed the damping performed in the computeAverage method.
		//			See the documentation at the beginning of the file for a motivation.
			void computeAverage(float a_damping) { avg = a_damping * avg + cur; }
			//void computeAverage(float a_damping) { avg = a_damping * (avg - cur) + cur; }

			void clear(void) { cur = 0; avg = 0; }
		};


		Data<uint32_t> entryCount;
		Data<tick_t> selfTicks;
		Data<tick_t> childTicks;


		tick_t totalTicksCur(void) const { return selfTicks.cur + childTicks.cur; }
		float totalTicksAvg(void) const { return selfTicks.avg + childTicks.avg; }

		void computeAverage(float a_damping) {
			entryCount.computeAverage(a_damping);
			selfTicks.computeAverage(a_damping);
			childTicks.computeAverage(a_damping);
		}

		void clearAll(void) {
			entryCount.clear();
			selfTicks.clear();
			childTicks.clear();
		}

		void clearCurrent(void) {
			entryCount.cur = 0;
			selfTicks.cur = 0;
			childTicks.cur = 0;
		}
	};


} // namespace Shiny

#endif

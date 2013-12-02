/*
 * memtracker.h
 *
 *  Created on: 08.05.2013
 *      Author: mrupp
 */

#ifndef MEMTRACKER_H_
#define MEMTRACKER_H_

#include "common/profiler/profiler.h"

namespace ug{

void UpdateTotalMem();
#ifdef UG_PROFILER_SHINY
size_t GetSelfMem(const Shiny::ProfileNode *p);
size_t GetTotalMem(const Shiny::ProfileNode *p);
#endif
bool EnableMemTracker(bool b);
bool IsMemTrackerEnabled();
void DisplayVacantMemory();
bool HasMemTracking();
}
#endif /* MEMTRACKER_H_ */

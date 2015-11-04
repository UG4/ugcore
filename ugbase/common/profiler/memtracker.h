/*
 * memtracker.h
 *
 *  Created on: 08.05.2013
 *      Author: mrupp
 *
 * You can enable the internal memtrack with
 * cmake -DPROFILE_MEMORY=ON ..
 *
 * If PROFILE_MEMORY=ON, global new and delete
 * operators are overwritten.
 * If INTERNAL_MEMTRACKER=OFF, global new and delete are untouched
 * and all other functions defined here lead to dummy functions.
 *
 * With EnableMemTracker(true) you start the memory profiling,
 * with EnableMemTracker(false) you stop it.
 * If you have a section you don't want to be measured, use this code
 * bool bMTEnabled = EnableMemTracker(false);
 * \code
 * CallToFunctionWhichShouldNotBeProfiled();
 * EnableMemTracker(bMTEnabled);
 * \endcode
 *
 * Memtracking can be output within
 * WriteProfileData("laplace.pdxml")
 * as part of the XML file and viewed with the ProfileViewer 1.32+
 * (http://gcsc.uni-frankfurt.de/Members/mrupp/shinyviewer/shinyprofileviewer)
 * It is also printed with console-commands like SetOutputProfileStats(true).
 *
 */

#ifndef MEMTRACKER_H_
#define MEMTRACKER_H_

#include "common/profiler/profiler.h"


namespace ug{

/**
 * The profiling process only calcs the selfMemory of each Profile Node.
 * This function calcs the total mem of each profile node.
 */
void UpdateTotalMem();

#ifdef UG_PROFILER_SHINY
size_t GetSelfMem(const Shiny::ProfileNode *p);
size_t GetTotalMem(const Shiny::ProfileNode *p);
#endif

/**
 * With EnableMemTracker(true) you start the memory profiling,
 * with EnableMemTracker(false) you stop it.
 * If you have a section you don't want to be measured, use this code
 * bool bMTEnabled = EnableMemTracker(false);
 * \code
 * CallToFunctionWhichShouldNotBeProfiled();
 * EnableMemTracker(bMTEnabled);
 * \endcode
 * EnableMemTracker(false) has to be called before pcl::Finalize();
 * @param b  if true, mem tracking is on, if false, it is off.
 * @return previous state.
 */
bool EnableMemTracker(bool b);

/**
 * @return true if currently we are memory tracking
 */
bool IsMemTrackerEnabled();


/**
 * @return true if -DINTERNAL_MEMTRACKER=ON, otherwise false
 * this function returns true if memory tracking is possible at all.
 */
bool HasMemTracking();


/**
 * Displays information about memory which has not been freed (yet).
 */
void DisplayVacantMemory();
}
#endif /* MEMTRACKER_H_ */

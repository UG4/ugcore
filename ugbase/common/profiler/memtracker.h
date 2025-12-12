/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/*
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


namespace ug {

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
#endif
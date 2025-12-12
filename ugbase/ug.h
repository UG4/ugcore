/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG__UG__
#define __H__UG__UG__

#include <string>

#include "common/profiler/profiler.h"
#include "common/ug_config.h"

/**
 * 	\brief the ug namespace
 *
 * Namespace for ug
 */
namespace ug {

////////////////////////////////////////////////////////////////////////
//	INFORMATION ON UG
///	Returns the version number of the current ug-version
/**	The string is formatted like this: "majorVersion.minorVersion.updateVersion"*/
UG_API std::string UGGetVersionString();


////////////////////////////////////////////////////////////////////////
//	INITIALISATION AND FINALISATION
///	initializes ug
/**	This method should be called at the beginning of main(...).
 *	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Init.
 *
 *	This method also sets the common paths in PathProvider.
 */
UG_API int UGInit(int *argcp, char ***argvp, int parallelOutputProcRank = -1);

///	Initializes the pathes of ug::PathProvider.
/**	Initializes the following pathes in ug::PathProvider:
 *	- ROOT_PATH
 *	- BIN_PATH
 *	- SCRIPT_PATH
 *	- DATA_PATH
 *	- GRID_PATH
 *	- PLUGIN_PATH
 *
 * Note: If you set a path before calling this method, it won't be overwritten.*/
UG_API bool InitPaths(const char* argv0);

///	Initializes the paths of ug::PathProvider.
/**	Initializes the following paths in ug::PathProvider relative to passed root path
 *	- ROOT_PATH
 *	- BIN_PATH
 *	- SCRIPT_PATH
 *	- DATA_PATH
 *	- GRID_PATH
 *	- PLUGIN_PATH
 */
UG_API void SetRootPath(const std::string& strRoot);
UG_API void SetRootPath(const char* c_strRoot);

///	Initializes the SCRIPT_PATH of ug::PathProvider.
UG_API void SetScriptPath(const std::string& strScript);
UG_API void SetScriptPath(const char* c_strScript);

///	Initializes the APPS_PATH of ug::PathProvider.
UG_API void SetAppsPath(const std::string& strApps);
UG_API void SetAppsPath(const char* c_strApps);

///	Initializes the PLUGIN_PATH of ug::PathProvider.
UG_API void SetPluginPath(const std::string& strPlugin);
UG_API void SetPluginPath(const char* c_strPlugin);

///	finalizes ug
/**	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Finalize.
 *
 *	Use ug::UGOutputProfileStatsOnExit to enable profiling output
 *	during finalize.
 */
UG_API int UGFinalize();

///	Calls UGFinalize and terminates the application.
/**	If the build-target is vrl and if no parallel build is performed, this method
 * throws an instance of SoftAbort*/
UG_API void UGForceExit();

///	Call with true, if profiling output is desired at the end of the show.
UG_API void UGOutputProfileStatsOnExit(bool bEnable);

///	Init (if UG_PLUGINS is set) embedded or non-shared plugins
UG_API bool UGInitPlugins();

///	sets a flag, that the current run shall be aborted during the next call of TerminateAbortedRun()
UG_API void AbortRun();
///	clears the abort-run-flag.
UG_API void ClearAbortRunFlag();
///	Terminates the current run if AbortRun() was called and the abort-run-flag is thus set to true.
UG_API void TerminateAbortedRun();

}//	end of namespace

#endif

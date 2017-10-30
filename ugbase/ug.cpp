/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
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
#include <stack>
#include <cstdlib>
#include <string>
#include "ug.h"
#include "common/error.h"
#include "common/log.h"
#include "common/util/path_provider.h"
#include "common/util/os_info.h"
#include "common/profiler/profiler.h"
#include "common/profiler/profile_node.h"

#include "common/profiler/memtracker.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif
#ifdef UG_BRIDGE
	#include "bridge/bridge.h"
#endif
#ifdef UG_PLUGINS
	#include "common/util/plugin_util.h"
#endif


/**	Current ug version */
//	ATTENTION: Please do not change the ug-version on your own!
//	If you changed something that requires a new version number, please contact
//	sreiter@gcsc.uni-frankfurt.de (for the moment...)
//	See docs/ug4/additional_pages/releases.doxygen for information on the different releases.
std::string gUGVersionString("4.0.2");


/** Tells whether profile-output is desired on exit.*/
//	only visible in this file!
static bool outputProfileStats = false;


namespace ug
{
	
static bool s_abortRun = false;


UG_API std::string UGGetVersionString()
{
	return gUGVersionString;
}


/**
 *  init app, script and data paths
 */
bool InitPaths(const char* argv0)
{
	PROFILE_FUNC();
	//The method currently only works if the path is explicitly specified
	//during startup or if UG4_ROOT is defined.

	//	extract the application path.
	char* ug4Root = getenv("UG4_ROOT");
	const char* pathSep = GetPathSeparator();

	std::string strRoot = "";

	if(ug4Root){
		strRoot = ug4Root;
	}
	else{
		std::string tPath = argv0;
		size_t pos = tPath.find_last_of(pathSep);
		
		if (pos != std::string::npos)
			tPath = tPath.substr(0, pos);
		else
			tPath = ".";

		strRoot = tPath + pathSep + "..";
	}

	if(!PathProvider::has_path(ROOT_PATH))
		PathProvider::set_path(ROOT_PATH, strRoot);
	if(!PathProvider::has_path(BIN_PATH))
		PathProvider::set_path(BIN_PATH, strRoot + pathSep + "bin");
	if(!PathProvider::has_path(SCRIPT_PATH))
		PathProvider::set_path(SCRIPT_PATH, strRoot + pathSep + "ugcore" + pathSep + "scripts");
	if(!PathProvider::has_path(PLUGIN_PATH))
		PathProvider::set_path(PLUGIN_PATH, strRoot + pathSep
										+ "bin" + pathSep + "plugins");
	if(!PathProvider::has_path(APPS_PATH))
		PathProvider::set_path(APPS_PATH, strRoot + pathSep + "apps");

//	log the paths
	UG_DLOG(MAIN, 1, "app path set to: " << PathProvider::get_path(BIN_PATH) <<
			std::endl << "script path set to: " << PathProvider::get_path(SCRIPT_PATH) <<
			std::endl);
/*
	if(!script::FileExists(PathProvider::get_path(BIN_PATH).c_str()) ||
	   !script::FileExists(PathProvider::get_path(SCRIPT_PATH).c_str()) ||
	   !script::FileExists(PathProvider::get_path(DATA_PATH).c_str()))
	{
		UG_LOG("WARNING: paths were not initialized correctly.\n");
		return false;
	}
*/
	return true;
}


/**
 *  init app, script and data paths for a given root path
 */
void SetRootPath(const std::string& strRoot)
{
	PROFILE_FUNC();
	const char* pathSep = GetPathSeparator();

	PathProvider::set_path(ROOT_PATH, strRoot);
	PathProvider::set_path(BIN_PATH, strRoot + pathSep + "bin");
	PathProvider::set_path(SCRIPT_PATH, strRoot + pathSep + "ugcore" + pathSep + "scripts");
	PathProvider::set_path(PLUGIN_PATH, strRoot + pathSep + "bin" + pathSep + "plugins");
	PathProvider::set_path(APPS_PATH, strRoot + pathSep + "apps");
}

////////////////////////////////////////////////////////////////////////
///	initializes ug
/**	This method should be called at the beginning of main(...).
 *	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Init.
 */
int UGInit(int *argcp, char ***argvp, int parallelOutputProcRank)
{
	PROFILE_FUNC();
	bool success = true;

	static bool firstCall = true;
	if (firstCall) {
		firstCall = false;

#ifdef UG_PARALLEL
		pcl::Init(argcp, argvp);
		GetLogAssistant().set_output_process(parallelOutputProcRank);
#endif

		success &= InitPaths((*argvp)[0]);

#ifdef UG_BRIDGE
		try{
			bridge::InitBridge();
		}
		catch(UGError& err)
		{
			success &= false;
			UG_LOG("ERROR in UGInit: InitBridge failed!\n");
		}
#endif

		if(UGInitPlugins() == false)
		{
			success &= false;
			UG_LOG("ERROR in UGInit: LoadPlugins failed!\n");
		}
	}

	// convert boolean success == true to int = 0.
	return !success;
}

int UGFinalizeNoPCLFinalize()
{
	EnableMemTracker(false);
	ug::GetLogAssistant().flush_error_log();
	
	if (outputProfileStats) {
		UG_LOG(std::endl);
	//	output the profiled data.
		PROFILER_UPDATE();

		if(GetLogAssistant().is_output_process()) {
			UG_LOG("\n");
#ifdef UG_PROFILER
			UG_LOG(ug::GetProfileNode(NULL)->call_tree());
#else
			PROFILER_OUTPUT();
#endif
		}

#ifdef UG_PROFILER
		//Shiny::ProfileManager::instance.destroy();
#endif
	}
	return 0;
}
////////////////////////////////////////////////////////////////////////
///	finalizes ug
/**	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Finalize.
 */
int UGFinalize()
{
	UGFinalizeNoPCLFinalize();

#ifdef UG_PARALLEL
	pcl::Finalize();
#endif

	return 0;
}

void UGForceExit()
{
	UG_LOG("--- ABORTING UG EXECUTION ---\n");

	#ifdef UG_PLUGINS
		// ? UnloadPlugins();
	#endif

	UGFinalizeNoPCLFinalize();

	#ifdef UG_PARALLEL
		if(pcl::NumProcs() > 1){
			// pcl::Abort will terminate execution
			pcl::Abort();
			// !!! this point is not reached !!!
		}
	#endif

	throw(SoftAbort("Exit forced by call to UGForceExit"));
}

void UGOutputProfileStatsOnExit(bool bEnable)
{
	outputProfileStats = bEnable;
}


#ifdef UG_PLUGINS
bool UGInitPlugins()
{
	return LoadPlugins(PathProvider::get_path(PLUGIN_PATH).c_str(), "ug4/", bridge::GetUGRegistry());
}
#else
bool UGInitPlugins()
{
	return true;
}
#endif


void AbortRun()
{
	s_abortRun = true;
}

void ClearAbortRunFlag()
{
	s_abortRun = false;
}

void TerminateAbortedRun()
{
	if(s_abortRun == true){
		s_abortRun = false;
		UGForceExit();
	}
}

} //end of namespace ug

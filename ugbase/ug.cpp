
#include <stack>
#include <cstdlib>
#include <string>
#include "ug.h"
#include "common/log.h"
#include "common/util/path_provider.h"
#include "common/util/os_info.h"
#include "common/profiler/profiler.h"
#include "common/profiler/profile_node.h"

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
		PathProvider::set_path(SCRIPT_PATH, strRoot + pathSep + "scripts");
	if(!PathProvider::has_path(DATA_PATH))
		PathProvider::set_path(DATA_PATH, strRoot + pathSep + "data");
	if(!PathProvider::has_path(GRID_PATH))
		PathProvider::set_path(GRID_PATH, strRoot + pathSep + "data" + pathSep + "grids");
	if(!PathProvider::has_path(PLUGIN_PATH))
		PathProvider::set_path(PLUGIN_PATH, strRoot + pathSep
										+ "bin" + pathSep + "plugins");
	if(!PathProvider::has_path(APPS_PATH))
		PathProvider::set_path(APPS_PATH, strRoot + pathSep + "apps");

//	log the paths
	UG_DLOG(MAIN, 1, "app path set to: " << PathProvider::get_path(BIN_PATH) <<
			std::endl << "script path set to: " << PathProvider::get_path(SCRIPT_PATH) <<
			std::endl << "data path set to: " << PathProvider::get_path(DATA_PATH) <<
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

////////////////////////////////////////////////////////////////////////
///	finalizes ug
/**	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Finalize.
 */
int UGFinalize()
{
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
	UGFinalize();
	exit(0);
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


} //end of namespace ug

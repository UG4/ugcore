/*
 * ug.cpp
 *
 *  Created on: 10.11.2010
 *      Author: marscher
 */

#include "ug.h"
#include "common/log.h"
#include "ug_script/ug_script.h"
#include "ug_bridge/ug_bridge.h"

#ifdef UG_PARALLEL
#include "pcl/pcl.h"
#endif

namespace ug {

//	we store the path in which scripts are located in this variable
//	we assume that scripts are located in ../scripts relative to ugshells path.

struct UGPaths {
	std::string APPS;
	std::string SCRIPTS;
	std::string DATA;
};

static struct UGPaths PATHS;

/**
 *  init app and data paths
 *  @return whether paths initialized correctly?
 */
static bool InitPaths(const char* argv0) {
	//TODO: on some systems argv0 does __not__ contain the absolute path to the process!
	//some ugly macros are needed.

	//	extract the application path.
	// UG_LOG("argv[0]: " << argv0 << endl);
	std::string tPath = argv0;
	size_t pos = tPath.find_last_of("/");
	if (pos == std::string::npos)
		pos = tPath.find_last_of("\\\\");
	if (pos != std::string::npos)
		PATHS.APPS = tPath.substr(0, pos);
	else
		PATHS.APPS = ".";

	PATHS.SCRIPTS = PATHS.APPS + "/../scripts";

	PATHS.DATA = std::string(PATHS.APPS);
	PATHS.DATA.append("/../data");

	UG_DLOG(MAIN, 0, "app path set to: " << PATHS.APPS <<
			std::endl << "script path set to: " << PATHS.SCRIPTS <<
			std::endl << "data path set to: " << PATHS.DATA << std::endl);

	if(!script::FileExists(PATHS.APPS.c_str()) ||
	   !script::FileExists(PATHS.SCRIPTS.c_str()) ||
	   !script::FileExists(PATHS.DATA.c_str())) {
		return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
///	initializes ug
/**	This method should be called at the beginning of main(...).
 *	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Init.
 */
//int UGInit(int argc, char* argv[], int parallelOutputProcRank)
int UGInit(int *argcp, char ***argvp, int parallelOutputProcRank) {
	//	make sure that things are only initialized once
	// todo: afaik static in d methods is per-cpp-file
	static bool firstCall = true;
	if (firstCall) {
		firstCall = false;
#ifdef UG_PARALLEL
//		pcl::Init(argc, argv);
		pcl::Init(argcp, argvp);
		pcl::SetOutputProcRank(parallelOutputProcRank);
#endif

		//	initialize ug-interfaces
		bridge::RegisterStandardInterfaces(bridge::GetUGRegistry());
	}

//	bool pathsCorrect = InitPaths(argv[0]);
	bool pathsCorrect = InitPaths((*argvp)[0]);
	if (!pathsCorrect)
		return -1;

	return 0;
}

/// returns the ug app path
const std::string& UGGetApplicationPath() {
	return PATHS.APPS;
}

/// returns the ug script path
const std::string& UGGetScriptPath()
{
	return PATHS.SCRIPTS;
}

/// returns the ug data path
const std::string& UGGetDataPath() {
	return PATHS.DATA;
}

////////////////////////////////////////////////////////////////////////
///	finalizes ug
/**	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Finalize.
 */
int UGFinalize(bool outputProfilerStats) {
	if (outputProfilerStats) {
	//	output the profiled data.
		PROFILER_UPDATE();
		#ifdef UG_PARALLEL
			if(pcl::IsOutputProc()) {
				UG_LOG("\n");
				PROFILER_OUTPUT();
			}
		#else
			UG_LOG("\n");
			PROFILER_OUTPUT();
		#endif
	}

#ifdef UG_PARALLEL
	pcl::Finalize();
#endif

	return 0;
}

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns its position in argv.
 *	If the parameter is not contained in argv, -1 is returned.
 */
int GetParamIndex(const char* param, int argc, char* argv[]) {
	for (int i = 0; i < argc; ++i) {
		if (strcmp(param, argv[i]) == 0) {
			return i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns true if it is found.
 */
bool FindParam(const char* param, int argc, char* argv[]) {
	return GetParamIndex(param, argc, argv) != -1;
}

////////////////////////////////////////////////////////////////////////
bool ParamToInt(int& iOut, const char* param, int argc, char* argv[]) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	iOut = atoi(argv[i + 1]);
	return true;
}

////////////////////////////////////////////////////////////////////////
bool ParamToString(char** strOut, const char* param, int argc, char* argv[]) {
	int i = GetParamIndex(param, argc, argv);
	if (i == -1 || i + 1 >= argc) {
		return false;
	}
	*strOut = argv[i + 1];
	return true;
}

} //end of namespace ug

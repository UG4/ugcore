// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m05 d31

#ifndef __H__UG__UG__
#define __H__UG__UG__
#include <string>

#include "common/profiler/profiler.h"

#include "common/common.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#include "lib_grid/lib_grid.h"
#include "node_tree/node_tree.h"
#include "ug_bridge/ug_bridge.h"
#include "ug_script/ug_script.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif


namespace ug
{

//	we store the path in which scripts are located in this variable
//	we assume that scripts are located in ../scripts relative to ugshells path.
static std::string APP_PATH;
static std::string SCRIPT_PATH;
static std::string DATA_PATH;

/**
 *  init app and data paths
 *  @return wether paths initialized correctly?
 */
inline static bool initPaths(const char* argv0) {
	//TODO: on some systems argv0 does __not__ contain the absolute path to the process!
	//some ugly macros are needed.


	//	extract the application path.
	cout << "argv[0]: " << argv0 << endl;
	string tPath = argv0;
	size_t pos = tPath.find_last_of("/");
	if(pos == string::npos) {
		pos = tPath.find_last_of("\\\\");
	}
	if(pos != string::npos) { // todo: this could be simply else?
		APP_PATH = tPath.substr(0, pos);
	}
	else {
		APP_PATH = ".";
	}

	SCRIPT_PATH = APP_PATH + "/../scripts";

	DATA_PATH=string(APP_PATH);
	DATA_PATH.append("/../data");

	LOG("app path set to: " << APP_PATH <<
			endl << "script path set to: " << SCRIPT_PATH <<
			endl << "data path set to: " << DATA_PATH << endl);

	if(!script::FileExists(APP_PATH.c_str()) || !script::FileExists(DATA_PATH.c_str()) || !script::FileExists(SCRIPT_PATH.c_str()) )  {
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
inline int UGInit(int argc, char* argv[], int parallelOutputProcRank = -1)
{
//	make sure that things are only initialized once
	// todo: afaik static in inlined methods is per-cpp-file
	static bool firstCall = true;
	if(firstCall){
		firstCall = false;
	#ifdef UG_PARALLEL
		pcl::Init(argc, argv);
		pcl::SetOutputProcRank(parallelOutputProcRank);
	#endif

	//	initialize ug-interfaces
		bridge::RegisterStandardInterfaces(bridge::GetUGRegistry());
	}

	bool pathsCorrect = initPaths(argv[0]);
	if(!pathsCorrect)
		return -1;

	return 0;
}

/// returns the ug app path
inline const std::string& UGGetApplicationPath() {
	return APP_PATH;
}
/// returns the ug script path
inline const std::string& UGGetScriptPath() {
	return SCRIPT_PATH;
}

/// returns the ug data path
inline const std::string& UGGetDataPath(){
	return DATA_PATH;
}

////////////////////////////////////////////////////////////////////////
///	finalizes ug
/**	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Finalize.
 */
inline int UGFinalize(bool outputProfilerStats = false)
{
	if(outputProfilerStats){
	//	output the profiled data.
		PROFILER_UPDATE();
		#ifdef UG_PARALLEL
			if(pcl::IsOutputProc()){
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
inline int GetParamIndex(const char* param, int argc, char* argv[])
{
	for(int i = 0; i < argc; ++i){
		if(strcmp(param, argv[i]) == 0){
			return i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns true if it is found.
 */
inline bool FindParam(const char* param, int argc, char* argv[])
{
	return GetParamIndex(param, argc, argv) != -1;
}

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter, and converts the
 *	associated value to an integer. Returns true if the parameter was
 *	found, false if not.
 */
inline bool ParamToInt(int& iOut, const char* param, int argc, char* argv[])
{
	int i = GetParamIndex(param, argc, argv);
	if(i==-1 || i+1 >= argc)
		return false;
	iOut = atoi(argv[i+1]);
	return true;
}



}//	end of namespace

#endif

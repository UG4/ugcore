// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m05 d31

#ifndef __H__UG__UG__
#define __H__UG__UG__

#include "common/profiler/profiler.h"

#include "common/common.h"
#include "lib_algebra/lib_algebra.h"
#include "lib_discretization/lib_discretization.h"
#include "lib_grid/lib_grid.h"
#include "node_tree/node_tree.h"
#include "ug_bridge/ug_bridge.h"
#include "ug_script/ug_script.h"
#include "bindings_vrl/bindings_vrl.h"

#ifdef UG_PARALLEL
	#include "pcl/pcl.h"
#endif


namespace ug
{
///	Returns ugs default registry.
inline bridge::Registry& GetUGRegistry()
{
	static bridge::Registry ugReg;
	return ugReg;
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
	static bool firstCall = true;
	if(firstCall){
		firstCall = false;
	#ifdef UG_PARALLEL
		pcl::Init(argc, argv);
		pcl::SetOutputProcRank(parallelOutputProcRank);
	#endif

	//	initialize ug-interfaces
		static bridge::Registry ugReg;
		bridge::RegisterStandardInterfaces(ugReg);
		script::SetScriptRegistry(&ugReg);
	}
	return 0;
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
/**	searches argv for the given parameter, and converts the
 *	associated value to an integer. Returns true if the parameter was
 *	found, false if not.
 */
inline bool ParamToInt(int& iOut, const char* param, int argc, char* argv[])
{
	for(int i = 0; i < argc; ++i){
		if(strcmp(param, argv[i]) == 0){
			if(argc > i + 1){
				iOut = atoi(argv[i+1]);
				return true;
			}
		}
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
/**	searches argv for the given parameter and returns true if it is found.
 */
inline bool FindParam(const char* param, int argc, char* argv[])
{
	for(int i = 0; i < argc; ++i){
		if(strcmp(param, argv[i]) == 0){
			return true;
		}
	}
	return false;
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


}//	end of namespace

#endif

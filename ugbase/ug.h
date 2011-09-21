/**
 * \file ug.h
 */
// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m05 d31

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
namespace ug
{

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

///	finalizes ug
/**	If ug has been compiled for parallel use (UG_PARALLEL is defined)
 *	then this method will internally call pcl::Finalize.
 *
 *	Use ug::UGOutputProfileStatsOnExit to enable profiling output
 *	during finalize.
 */
UG_API int UGFinalize();

///	Calls UGFinalize and terminates the application.
UG_API void UGForceExit();

///	Call with true, if profiling output is desired at the end of the show.
UG_API void UGOutputProfileStatsOnExit(bool bEnable);

}//	end of namespace

#endif

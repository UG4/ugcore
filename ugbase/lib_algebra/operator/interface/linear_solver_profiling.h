
#ifndef LINEAR_SOLVER_PROFILING_H_
#define LINEAR_SOLVER_PROFILING_H_

#include "common/profiler/profiler.h"
#define PROFILE_LS
#ifdef PROFILE_LS
	#define LS_PROFILE_FUNC()		PROFILE_FUNC_GROUP("LinearSolver algebra")
	#define LS_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "LinearSolver algebra")
	#define LS_PROFILE_END(name)		PROFILE_END_(name)
#else
	#define LS_PROFILE_FUNC()
	#define LS_PROFILE_BEGIN(name)
	#define LS_PROFILE_END(name)
#endif

#endif /* LINEAR_SOLVER_PROFILING_H_ */

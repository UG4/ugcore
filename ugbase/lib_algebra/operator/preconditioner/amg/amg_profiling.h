/*
 * amg_profiling.h
 *
 *  Created on: 11.04.2011
 *      Author: mrupp
 */

#ifndef AMG_PROFILING_H_
#define AMG_PROFILING_H_


#define PROFILE_AMG
#ifdef PROFILE_AMG
	#include  "common/profiler/profiler.h"
	#define AMG_PROFILE_FUNC()			PROFILE_FUNC()
	#define AMG_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define AMG_PROFILE_END()			PROFILE_END()
#else
	#define AMG_PROFILE_FUNC()
	#define AMG_PROFILE_BEGIN(name)
	#define AMG_PROFILE_END()
#endif

#endif /* AMG_PROFILING_H_ */

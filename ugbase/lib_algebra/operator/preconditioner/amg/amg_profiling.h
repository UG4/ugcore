/*
 * amg_profiling.h
 *
 *  Created on: 11.04.2011
 *      Author: mrupp
 */

#ifndef __H__LIB_ALGEBRA__AMG_PROFILING_H_
#define __H__LIB_ALGEBRA__AMG_PROFILING_H_


#define PROFILE_AMG
#ifdef PROFILE_AMG
	#include  "common/profiler/profiler.h"
	#define AMG_PROFILE_FUNC()			PROFILE_FUNC()
	#define AMG_PROFILE_BEGIN(name)		PROFILE_BEGIN(name)
	#define AMG_PROFILE_NEXT(name)		PROFILE_END(); PROFILE_BEGIN(name)
	#define AMG_PROFILE_END()			PROFILE_END()
	#define AMG_BPROFILE_BEGIN(name, string) PROFILE_BEGIN(name); stopwatch SW; const char *profileString = string; SW.start();
	#define AMG_BPROFILE_END() PROFILE_END(); UG_DLOG(LIB_ALG_AMG, 1, string << " took " << SW.ms() << " ms\n");
	#define AMG_BPROFILE_NEXT(name, string) AMG_BPROFILE_END(); PROFILE_BEGIN(name); SW.start();

#else
	#define AMG_PROFILE_FUNC()
	#define AMG_PROFILE_BEGIN(name)
	#define AMG_PROFILE_NEXT(name)
	#define AMG_PROFILE_END()
	#define AMG_BPROFILE_BEGIN(name, string)
	#define AMG_BPROFILE_NEXT(name, string)
	#define AMG_BPROFILE_END()
#endif

#endif /* __H__LIB_ALGEBRA__AMG_PROFILING_H_ */

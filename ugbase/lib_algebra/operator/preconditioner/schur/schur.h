/*
 * schur.h
 *
 *  Created on: 18.12.2013
 *      Author: anaegel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__


#ifdef UG_PARALLEL

#include "common/profiler/profiler.h"
#include "common/log.h"

#define PROFILE_SCHUR
#ifdef PROFILE_SCHUR
	#define SCHUR_PROFILE_FUNC()			PROFILE_FUNC_GROUP("algebra schur")
	#define SCHUR_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "algebra schur")
	#define SCHUR_PROFILE_END_(name)			PROFILE_END_(name)
#else
	#define SCHUR_PROFILE_FUNC()
	#define SCHUR_PROFILE_BEGIN(name)
	#define SCHUR_PROFILE_END_(name)
#endif

namespace ug{

extern DebugID SchurDebug;

} // end namespace ug

#include "schur_complement_operator.h"
#include "schur_precond.h"

#endif /* UG_PARALLEL */

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__SCHUR__ */

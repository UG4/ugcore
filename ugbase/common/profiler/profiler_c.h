#ifndef __H__UG__COMMON__PROFILER_C__
#define __H__UG__COMMON__PROFILER_C__

#ifdef __H__UG__COMMON__PROFILER__
#error "can't use profiler_c.h and profiler.h"
#endif

#ifdef UG_PROFILER


#ifdef UG_PROFILER_SHINY

//	this is just a wrapper-include for the shiny-profiler by Aidin Abedi
	#define SHINY_PROFILER TRUE
	#include "src/ShinyManager.h"
	#include "src/ShinyNode.h"


	#define C_PROFILE_BEGIN_PARAMS(id, name, group, file, line)			\
															\
		static Shiny::ProfileZone __ShinyZone_##id = {		\
			NULL, Shiny::ProfileZone::STATE_HIDDEN, name, \
			group, file, line,	\
			{ { 0, 0 }, { 0, 0 }, { 0, 0 } }				\
		};													\
		{													\
			static Shiny::ProfileNodeCache cache =			\
				&Shiny::ProfileNode::_dummy;				\
															\
			Shiny::ProfileManager::instance._beginNode(&cache, &__ShinyZone_##id);\
		}


	#define C_PROFILE_BEGIN_GROUP(name, group)						\
		C_PROFILE_BEGIN_PARAMS(apn_##name, #name, group, __FILE__, __LINE__)

	#define C_PROFILE_BEGIN(name)						\
		C_PROFILE_BEGIN_GROUP(name, NULL)

	#define C_PROFILE_END()														\
		Shiny::ProfileManager::instance._endCurNode()

	#define C_PROFILE_FUNC_BEGIN() \
		C_PROFILE_BEGIN(__FUNCTION__)
	#define PROFILE_FUNC_GROUP_BEGIN(groups) \
		C_PROFILE_BEGIN_GROUP(__FUNCTION__, groups)

#else
#error "not defined for C"
#endif // UG_PROFILER_SHINY

#else
	#define C_PROFILE_BEGIN(name)
	#define C_PROFILE_BEGIN_GROUP(name, groups)
	#define C_PROFILE_END()
	#define C_PROFILE_FUNC_BEGIN()
	#define C_PROFILE_FUNC_GROUP_BEGIN(groups)

#endif // UG_PROFILER

#endif	// __H__UG__COMMON__PROFILER__

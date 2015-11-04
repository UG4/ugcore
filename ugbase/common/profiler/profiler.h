#ifndef __H__UG__COMMON__PROFILER__
#define __H__UG__COMMON__PROFILER__

// if cpu-frequency adaption is enabled
#ifdef UG_CPU_FREQ
	#include "freq_adapt.h"
	#define CPU_FREQ_BEGIN_AUTO_END(id, file, line)    \
			static unsigned long __freq__##id = FreqAdaptValues::freq(file, line); \
			AutoFreqAdaptNode __node__freq__##id((__freq__##id));

	#define CPU_FREQ_END()  \
			FreqAdaptNodeManager::release_latest();
#else
	#define CPU_FREQ_BEGIN_AUTO_END(name, file, line)
	#define CPU_FREQ_END()
#endif


// if some profiler is enabled
#ifdef UG_PROFILER

#include <vector>
#include "profilenode_management.h"
#include "shiny_call_logging.h"


#ifdef UG_PROFILER_SHINY

//	this is just a wrapper-include for the shiny-profiler by Aidin Abedi
	#define SHINY_PROFILER TRUE
	#include "src/ShinyManager.h"
	#include "src/ShinyNode.h"


	/**	Helper makro used in PROFILE_BEGIN and PROFILE_FUNC.*/
	#define PROFILE_BEGIN_AUTO_END(id, name, group, file, line)			\
															\
		CPU_FREQ_BEGIN_AUTO_END(id, file, line); 			\
		AutoProfileNode	id;									\
		static Shiny::ProfileZone __ShinyZone_##id = {		\
			NULL, Shiny::ProfileZone::STATE_HIDDEN, name, 	\
			group, file, line,								\
			{ { 0, 0 }, { 0, 0 }, { 0, 0 } }				\
		};													\
		{													\
			static Shiny::ProfileNodeCache cache =			\
				&Shiny::ProfileNode::_dummy;				\
															\
			Shiny::ProfileManager::instance._beginNode(&cache, &__ShinyZone_##id);\
		}\
		PROFILE_LOG_CALL_START()


	/**	Creates a new profile-environment with the given name.
	 * Note that the profiled section automatically ends when the current
	 * ends.
	 */
	#define PROFILE_BEGIN(name)						\
			PROFILE_BEGIN_AUTO_END(apn_##name, #name, NULL, __FILE__, __LINE__)

	/**	Ends profiling of the latest PROFILE_BEGIN section.*/
	#define PROFILE_END()							\
			ProfileNodeManager::release_latest(); \
			CPU_FREQ_END();

	/**	Profiles the whole function*/
	#define PROFILE_FUNC()										\
			PROFILE_BEGIN_AUTO_END(__ShinyFunction, __FUNCTION__, NULL, __FILE__, __LINE__)

	#define PROFILE_BEGIN_GROUP(name, group)					\
		PROFILE_BEGIN_AUTO_END(apn_##name, #name, group, __FILE__, __LINE__)

	#define PROFILE_FUNC_GROUP(group)										\
			PROFILE_BEGIN_AUTO_END(__ShinyFunction, __FUNCTION__, group, __FILE__, __LINE__)

	/**	Performs update on the profiler (call before output)*/
	#define PROFILER_UPDATE									\
		Shiny::ProfileManager::instance.update

	/**	Outputs the profile-times*/
	#define PROFILER_OUTPUT									\
		Shiny::ProfileManager::instance.output

#endif // UG_PROFILER_SHINY


#ifdef UG_PROFILER_SCALASCA
	#include "epik_user.h"
	#include <ostream>

	#define PROFILE_STRINGIFY(x) #x
	#define PROFILE_TOSTRING(x) PROFILE_STRINGIFY(x)

	/**	Creates a new profile-environment with the given name.
	 * Note that the profiled section automatically ends when the current ends.
	 */
	#define PROFILE_BEGIN(name)	\
		EPIK_USER_REG(__##name, PROFILE_TOSTRING(name));	\
		EPIK_USER_START(__##name);	\
		AutoProfileNode	apn_##name(__##name);								\

	/**	Ends profiling of the latest PROFILE_BEGIN section.*/
	#define PROFILE_END()										\
			ProfileNodeManager::release_latest()

	/**	Profiles the whole function*/
	#define PROFILE_FUNC()										\
			EPIK_TRACER(__FUNCTION__)

	#define PROFILE_BEGIN_GROUP(name, group)					\
			PROFILE_BEGIN(name)

	#define PROFILE_FUNC_GROUP(group)							\
			PROFILE_FUNC()

	namespace ProfilerDummy{
		inline void Update(float a = 0.0f)			{}
		inline bool Output(const char *a = NULL)	{return false;}
		inline bool Output(std::ostream &a)			{return false;}
	}

	#define PROFILER_UPDATE	ProfilerDummy::Update
	#define PROFILER_OUTPUT	ProfilerDummy::Output

#endif // UG_PROFILER_SCALASCA

#ifdef UG_PROFILER_VAMPIR
	#include "vt_user.h"
	#include <ostream>

	#define PROFILE_STRINGIFY(x) #x
	#define PROFILE_TOSTRING(x) PROFILE_STRINGIFY(x)

	/**	Creates a new profile-environment with the given name.
	 * Note that the profiled section automatically ends when the current ends.
	 */
	#define PROFILE_BEGIN(name)	\
			VT_USER_START(PROFILE_TOSTRING(name));	\
			AutoProfileNode	apn_##name(PROFILE_TOSTRING(name));	\

	/**	Ends profiling of the latest PROFILE_BEGIN section.*/
	#define PROFILE_END()										\
			ProfileNodeManager::release_latest()

	/**	Profiles the whole function*/
	#define PROFILE_FUNC()										\
			VT_TRACER((char*)__FUNCTION__)

	#define PROFILE_BEGIN_GROUP(name, group)					\
			PROFILE_BEGIN(name)

	#define PROFILE_FUNC_GROUP(group)							\
			PROFILE_FUNC()

	namespace ProfilerDummy{
		inline void Update(float a = 0.0f)			{}
		inline bool Output(const char *a = NULL)	{return false;}
		inline bool Output(std::ostream &a)			{return false;}
	}

	#define PROFILER_UPDATE	ProfilerDummy::Update
	#define PROFILER_OUTPUT	ProfilerDummy::Output

#endif // UG_PROFILER_VAMPIR

#ifdef UG_PROFILER_SCOREP
	#include <scorep/SCOREP_User.h>
	#include <ostream>

	#define PROFILE_STRINGIFY(x) #x
	#define PROFILE_TOSTRING(x) PROFILE_STRINGIFY(x)

	/**	Creates a new profile-environment with the given name.
	 * Note that the profiled section automatically ends when the current ends.
	 */
	#define PROFILE_BEGIN(name)	\
			SCOREP_USER_REGION_DEFINE( __scorephandle__##name )								\
			SCOREP_USER_REGION_BEGIN( __scorephandle__##name, PROFILE_TOSTRING(name),			\
			                          SCOREP_USER_REGION_TYPE_COMMON ) 			\
			AutoProfileNode	apn_##name(__scorephandle__##name);

	/**	Ends profiling of the latest PROFILE_BEGIN section.*/
	#define PROFILE_END()										\
			ProfileNodeManager::release_latest()

	/**	Profiles the whole function*/
	#define PROFILE_FUNC()										\
			SCOREP_USER_REGION(__FUNCTION__, SCOREP_USER_REGION_TYPE_FUNCTION )

	#define PROFILE_BEGIN_GROUP(name, group)					\
			PROFILE_BEGIN(name)

	#define PROFILE_FUNC_GROUP(group)							\
			PROFILE_FUNC()

	namespace ProfilerDummy{
		inline void Update(float a = 0.0f)			{}
		inline bool Output(const char *a = NULL)	{return false;}
		inline bool Output(std::ostream &a)			{return false;}
	}

	#define PROFILER_UPDATE	ProfilerDummy::Update
	#define PROFILER_OUTPUT	ProfilerDummy::Output

#endif // UG_PROFILER_SCOREP

#define PROFILE_END_(name) \
			assert(&(apn_##name) == ProfileNodeManager::inst().m_nodes.top());	\
			struct apn_already_ended_##name { } ; \
			PROFILE_END();

#else
	#include <ostream>

	namespace ProfilerDummy{
		inline void Update(float a = 0.0f)			{}
		inline bool Output(const char *a = NULL)	{return false;}
		inline bool Output(std::ostream &a)			{return false;}
	}

//	Empty macros if UG_PROFILER == false
	#define PROFILE_BEGIN(name) 				CPU_FREQ_BEGIN_AUTO_END(apn_##name, __FILE__, __LINE__);
	#define PROFILE_BEGIN_GROUP(name, groups) 	PROFILE_BEGIN(name)
	#define PROFILE_END() 						CPU_FREQ_END();
	#define PROFILE_FUNC() 						CPU_FREQ_BEGIN_AUTO_END(apn##__FUNCTION__, __FILE__, __LINE__);
	#define PROFILE_FUNC_GROUP(groups) 			PROFILE_FUNC()
	#define PROFILER_UPDATE	ProfilerDummy::Update
	#define PROFILER_OUTPUT	ProfilerDummy::Output

	#define PROFILE_END_(name)

#endif // UG_PROFILER

#endif	// __H__UG__COMMON__PROFILER__

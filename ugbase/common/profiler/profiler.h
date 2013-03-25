//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d21

#ifndef __H__UG__COMMON__PROFILER__
#define __H__UG__COMMON__PROFILER__

#ifdef UG_PROFILER

#include <stack>

#ifdef UG_PROFILER_SCOREP
	#include <scorep/SCOREP_User.h>
#endif

class AutoProfileNode;

class ProfileNodeManager
{
	public:
		static void add(AutoProfileNode* node);
		static void release_latest();

	private:
		ProfileNodeManager();
		~ProfileNodeManager();
		static ProfileNodeManager& inst();

	private:
		std::stack<AutoProfileNode*>	m_nodes;
};

class AutoProfileNode
{
	friend class ProfileNodeManager;

	public:
#ifdef UG_PROFILER_SHINY
		AutoProfileNode();
#endif
#if defined(UG_PROFILER_SCALASCA) || defined(UG_PROFILER_VAMPIR)
		AutoProfileNode(const char* name);
#endif
#ifdef UG_PROFILER_SCOREP
		AutoProfileNode(SCOREP_User_RegionHandle name);
#endif
		~AutoProfileNode();

	private:
		void release();
		inline bool is_active()		{return m_bActive;}

	private:
		bool m_bActive;
#if defined(UG_PROFILER_SCALASCA) || defined(UG_PROFILER_VAMPIR)
		const char* m_pName;
#endif
#ifdef UG_PROFILER_SCOREP
		SCOREP_User_RegionHandle m_pHandle;
#endif
};

#ifdef UG_PROFILER_SHINY

//	this is just a wrapper-include for the shiny-profiler by Aidin Abedi
	#define SHINY_PROFILER TRUE
	#include "src/ShinyManager.h"

	/**	Helper makro used in PROFILE_BEGIN and PROFILE_FUNC.*/
	#define PROFILE_BEGIN_AUTO_END(id, name, group, file, line)			\
															\
		AutoProfileNode	id;								\
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

	/**	Creates a new profile-environment with the given name.
	 * Note that the profiled section automatically ends when the current
	 * ends.
	 */
	#define PROFILE_BEGIN(name)						\
			PROFILE_BEGIN_AUTO_END(apn_##name, #name, NULL, __FILE__, __LINE__)

	/**	Ends profiling of the latest PROFILE_BEGIN section.*/
	#define PROFILE_END()							\
			ProfileNodeManager::release_latest()

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
			AutoProfileNode	apn_##name(name);

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

#else
	#include <ostream>

	namespace ProfilerDummy{
		inline void Update(float a = 0.0f)			{}
		inline bool Output(const char *a = NULL)	{return false;}
		inline bool Output(std::ostream &a)			{return false;}
	}

//	Empty macros if UG_PROFILER == false
	#define PROFILE_BEGIN(name)
	#define PROFILE_BEGIN_GROUP(name, groups)
	#define PROFILE_END()
	#define PROFILE_FUNC()
	#define PROFILE_FUNC_GROUP(groups)
	#define PROFILER_UPDATE	ProfilerDummy::Update
	#define PROFILER_OUTPUT	ProfilerDummy::Output

#endif // UG_PROFILER

#endif	// __H__UG__COMMON__PROFILER__

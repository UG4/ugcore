//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m10 d21

#ifndef __H__UG__COMMON__PROFILER__
#define __H__UG__COMMON__PROFILER__

//	this is just a wrapper-include for the shiny-profiler by Aidin Abedi

//	To enable or disable the profiler the following define can be set.
//	It is preferable to do this via cmake options or similar.
//#define UG_PROFILER
#ifdef UG_PROFILER
	#define SHINY_PROFILER TRUE
#else
	#define SHINY_PROFILER FALSE
#endif

#ifdef UG_PROFILER
	#include <stack>
	#include "src/ShinyManager.h"

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
			AutoProfileNode();
			~AutoProfileNode();

		private:
			inline void deactivate()	{m_bActive = false;}

		private:
			bool m_bActive;

	};


	/**	Helper makro used in PROFILE_BEGIN and PROFILE_FUNC.*/
	#define PROFILE_BEGIN_HELPER(id, name, group, file, line)			\
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
			PROFILE_BEGIN_HELPER(apn_##name, #name, NULL, __FILE__, __LINE__)

	/**	Ends profiling of the latest PROFILE_BEGIN section.*/
	#define PROFILE_END()							\
			ProfileNodeManager::release_latest()

	/**	Profiles the whole function*/
	#define PROFILE_FUNC()										\
			PROFILE_BEGIN_HELPER(__ShinyFunction, __FUNCTION__, NULL, __FILE__, __LINE__)

	#define PROFILE_BEGIN_GROUP(name, group)					\
		PROFILE_BEGIN_HELPER(apn_##name, #name, group, __FILE__, __LINE__)

	#define PROFILE_FUNC_GROUP(group)										\
			PROFILE_BEGIN_HELPER(__ShinyFunction, __FUNCTION__, group, __FILE__, __LINE__)

	/**	Performs update on the profiler (call before output)*/
	#define PROFILER_UPDATE									\
		Shiny::ProfileManager::instance.update

	/**	Outputs the profile-times*/
	#define PROFILER_OUTPUT									\
		Shiny::ProfileManager::instance.output

#else
	#include <iostream>

	namespace ShinyDummy{
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
	#define PROFILER_UPDATE	ShinyDummy::Update
	#define PROFILER_OUTPUT	ShinyDummy::Output

#endif // UG_PROFILER

#endif	// __H__UG__COMMON__PROFILER__

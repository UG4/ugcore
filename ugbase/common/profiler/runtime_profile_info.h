// Lua Profiling (mrupp)

#ifndef __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__
#define __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__

#include "profiler.h"

#ifdef UG_PROFILER

/**
 * Class storing Profile information, only known at runtime (e.g. strings build
 * at runtime in lua binding). This class can be used to also profile such
 * parts of the application. To do so, this class stores a string pointer to name,
 * group, file etc. that remains valid until the end of the class (and often
 * program). Therefore, it can be used like the static const char* identifiers
 * that are commonly used for profiling.
 *
 * RuntimeProfileInfo only copies the name and the group if you tell it
 * to! So pay attention if pointers you hand over to RuntimeProfileInfo
 * are persistent.
 */
class RuntimeProfileInfo
{
	public:
		RuntimeProfileInfo(const char* name = NULL, bool bCopyName = false,
		                          const char* groups = NULL, bool bCopyGroup = false,
		                          const char* file = NULL, bool bCopyFile = false,
		                          int line = 0);

		~RuntimeProfileInfo();

		inline void beginNode()
		{
			Shiny::ProfileManager::instance._beginNode(&profilerCache, &profileInformation);
		}

		inline void endNode()
		{
			Shiny::ProfileManager::instance._endCurNode();
		}

	private:
		Shiny::ProfileZone profileInformation;
		Shiny::ProfileNodeCache profilerCache;
		bool bNameCopied;
		bool bGroupCopied;
		bool bFileCopied;
		const char* pName;
		const char* pGroup;
		const char* pFile;
		int iLine;
};

typedef RuntimeProfileInfo* pRuntimeProfileInfo;
#endif


#endif /* __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__ */

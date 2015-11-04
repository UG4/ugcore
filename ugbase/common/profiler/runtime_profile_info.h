
#ifndef __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__
#define __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__

#ifdef UG_PROFILER

#include "profiler.h"
#include "common/log.h"
#include "common/util/stringify.h"

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
#ifdef UG_PROFILER_SHINY
		Shiny::ProfileManager::instance._beginNode(&profilerCache, &profileInformation);
		PROFILE_LOG_CALL_START();
#endif
#ifdef UG_PROFILER_SCALASCA
		EPIK_USER_START(pName);
#endif
#ifdef UG_PROFILER_VAMPIR
		VT_USER_START((char*)pName);
#endif
#ifdef UG_PROFILER_SCOREP
		SCOREP_USER_REGION_BEGIN( m_pHandle, pName,
								  SCOREP_USER_REGION_TYPE_COMMON )
#endif
	}

	inline void endNode()
	{
#ifdef UG_PROFILER_SHINY
		Shiny::ProfileManager::instance._endCurNode();
		PROFILE_LOG_CALL_END();
#endif
#ifdef UG_PROFILER_SCALASCA
		EPIK_USER_END(pName);
#endif
#ifdef UG_PROFILER_VAMPIR
		VT_USER_END((char*)pName);
#endif
#ifdef UG_PROFILER_SCOREP
		SCOREP_USER_REGION_END(m_pHandle);
#endif
	}


	const char *name() const
	{
		return pName;
	}

	const char *group() const
	{
		return pGroup;
	}

	const char *file() const
	{
		return pFile;
	}

	int line() const
	{
		return iLine;
	}

	private:
		bool bNameCopied;
		bool bGroupCopied;
		bool bFileCopied;
		const char* pName;
		const char* pGroup;
		const char* pFile;
		int iLine;

#ifdef UG_PROFILER_SHINY
		Shiny::ProfileZone profileInformation;
		Shiny::ProfileNodeCache profilerCache;
#endif
#ifdef UG_PROFILER_SCOREP
		SCOREP_User_RegionHandle m_pHandle;
#endif
};

static inline std::ostream& operator << (std::ostream& os, const RuntimeProfileInfo &pi)
{
	os << "RuntimeProfileInfo name=" << pi.name() << " group=" << pi.group() << " @ " << pi.file() << ":" << pi.line();
	return os;
}

typedef RuntimeProfileInfo* pRuntimeProfileInfo;
#endif


#endif /* __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__ */

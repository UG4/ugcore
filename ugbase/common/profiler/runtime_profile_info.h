/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */


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
	RuntimeProfileInfo(const char* name = nullptr, bool bCopyName = false,
					   const char* groups = nullptr, bool bCopyGroup = false,
					   const char* file = nullptr, bool bCopyFile = false,
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

using pRuntimeProfileInfo = RuntimeProfileInfo*;
#endif


#endif /* __H__UG__COMMON__PROFILER__RUNTIME_PROFILE_INFO__ */

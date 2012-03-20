#ifndef __H__UG__COMMON__DYNAMIC_PROFILING__
#define __H__UG__COMMON__DYNAMIC_PROFILING__
#include "profiler.h"

#ifdef UG_PROFILER
// Lua Profiling (mrupp)

struct DynamicProfileInformation
{
	DynamicProfileInformation(const char*name=NULL, bool bCopy=false)
	{
		Shiny::ProfileZone pi = {NULL, Shiny::ProfileZone::STATE_HIDDEN, NULL, { { 0, 0 }, { 0, 0 }, { 0, 0 } }};
		profileInformation = pi;
		profilerCache =	&Shiny::ProfileNode::_dummy;
		bCopied = false;
		if(name)
		{
			if(bCopy)
				set_name(name);
			else
				set_name_and_copy(name);
		}
	}
	~DynamicProfileInformation()
	{
		if(bCopied && profileInformation.name) delete[] profileInformation.name;
	}
	Shiny::ProfileZone profileInformation;
	Shiny::ProfileNodeCache profilerCache;
	bool bCopied;

	bool is_initialised()
	{
		return profileInformation.name != NULL;
	}

	void set_name_and_copy(const char*name)
	{
		if(bCopied && profileInformation.name) delete[] profileInformation.name;
		char *p = new char[strlen(name)+1];
		strcpy(p, name);
		name = p;
		bCopied = true;
		profileInformation.name = name;
	}

	inline void set_name(const char*name)
	{
		profileInformation.name = name;
	}

	inline void beginNode()
	{
		Shiny::ProfileManager::instance._beginNode(&profilerCache, &profileInformation);
	}

	static inline void endCurNode()
	{
		Shiny::ProfileManager::instance._endCurNode();
	}

	inline bool isCurNode()
	{
		return Shiny::ProfileManager::instance._curNode == profilerCache;
	}
};

typedef DynamicProfileInformation * pDynamicProfileInformation;
#endif

#endif

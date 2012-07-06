#ifndef __H__UG__COMMON__DYNAMIC_PROFILING__
#define __H__UG__COMMON__DYNAMIC_PROFILING__
#include "profiler.h"

#ifdef UG_PROFILER
// Lua Profiling (mrupp)

/**
 * Use this version to ProfileFunction of which you don't know the name at compile time.
 * ex:
 * ... some function
 * {
 * 		static DynamicProfileInformation dpi;
 * 		if(!dpi.is_initialized()) { a.init(myname, true, "mygroup", false); }
 * 		dpi.beginNode();
 * 		....
 * 		dpi.endCurNode();
 * 	}
 *
 * DynamicProfileInformation only copies the name and the group if you tell it to!
 * So pay attention if pointers you hand over to DynamicProfileInformation are persitent.
 */
struct DynamicProfileInformation
{
	DynamicProfileInformation(const char*name=NULL, bool bCopy=false, const char *groups=NULL, bool bCopyGroup=false)
	{
		Shiny::ProfileZone pi = {NULL, Shiny::ProfileZone::STATE_HIDDEN, NULL, NULL,
				{ { 0, 0 }, { 0, 0 }, { 0, 0 } }};
		profileInformation = pi;
		profilerCache =	&Shiny::ProfileNode::_dummy;
		bCopied = false;
		bgGroupCopied = false;
		init(name, bCopyName, group, bCopyGroup);
	}
	~DynamicProfileInformation()
	{
		if(bCopied && profileInformation.name) delete[] profileInformation.name;
	}

	bool is_initialised()
	{
		return profileInformation.name != NULL;
	}

	void init(const char*name=NULL, bool bCopy=false, const char*group=NULL, bool bCopyGroup=false)
	{
		if(bCopied && profileInformation.name) delete[] profileInformation.name;
		if(bGroupCopied && profileInformation.group) delete[] profileInformation.group;
		if(bCopy)
		{
			char *p= new char[strlen(name)+1];
			strcpy(p, name)
			bCopied = true;
			name = p;
		}
		if(bCopyGroup)
		{
			char *p = new char[strlen(group)+1];
			strcpy(p, group)
			bGroupCopied = true;
			group = p;
		}

		profileInformation.name = name;
		profileInformation.group = group;
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
private:
	Shiny::ProfileZone profileInformation;
	Shiny::ProfileNodeCache profilerCache;
	bool bCopied;
	bool bGroupCopied;
};

typedef DynamicProfileInformation * pDynamicProfileInformation;
#endif


#endif

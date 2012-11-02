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
	DynamicProfileInformation(
			const char*name=NULL, bool bCopyName=false, 
			const char *groups=NULL, bool bCopyGroup=false,
			const char *file=NULL, bool bCopyFile=false, 
			int line=0)
	{
		Shiny::ProfileZone pi = {NULL, Shiny::ProfileZone::STATE_HIDDEN, NULL, NULL, NULL, 0,
				{ { 0, 0 }, { 0, 0 }, { 0, 0 } }};
		profileInformation = pi;
		profilerCache =	&Shiny::ProfileNode::_dummy;
		bCopied = false;
		bGroupCopied = false;
		init(name, bCopyName, groups, bCopyGroup, file, bCopyFile, line);
	}
	~DynamicProfileInformation()
	{
		if(bCopied && profileInformation.name) delete[] profileInformation.name;
	}

	bool is_initialised()
	{
		return profileInformation.name != NULL;
	}

	void init(const char*name=NULL, bool bCopyName=false, const char*groups=NULL, bool bCopyGroup=false,
	const char *file=NULL, bool bCopyFile=false, int line=0)
	{
		if(bCopied && profileInformation.name) delete[] profileInformation.name;
		if(bGroupCopied && profileInformation.groups) delete[] profileInformation.groups;
		if(bCopyName)
		{
			char *p= new char[strlen(name)+1];
			strcpy(p, name);
			bCopied = true;
			name = p;
		}
		if(bCopyGroup)
		{
			char *p = new char[strlen(groups)+1];
			strcpy(p, groups);
			bGroupCopied = true;
			groups = p;
		}
		if(bCopyFile)
		{
			char *p = new char[strlen(file)+1];
			strcpy(p, file);
			bGroupCopied = true;
			file = p;
		}

		profileInformation.name = name;
		profileInformation.groups = groups;
		profileInformation.file = file;
		profileInformation.line = line;
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

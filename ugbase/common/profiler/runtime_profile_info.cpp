
#include "runtime_profile_info.h"
#include "string.h"

RuntimeProfileInfo::RuntimeProfileInfo(
		const char* name, bool bCopyName,
		const char* groups, bool bCopyGroup,
		const char* file, bool bCopyFile,
		int line)
: bNameCopied(false), bGroupCopied(false), bFileCopied(false),
  pName(NULL), pGroup(NULL), pFile(NULL)
{
	if(bCopyName) {
		char* p = new char[strlen(name)+1];
		strcpy(p, name);
		pName = p;
		bNameCopied = true;
	} else {
		pName = name;
	}

	if(bCopyGroup)
	{
		char* p = new char[strlen(groups)+1];
		strcpy(p, groups);
		pGroup = p;
		bGroupCopied = true;
	} else {
		pGroup = groups;
	}

	if(bCopyFile)
	{
		char* p = new char[strlen(file)+1];
		strcpy(p, file);
		pFile = p;
		bFileCopied = true;
	} else {
		pFile = file;
	}
	iLine = line;

#ifdef UG_PROFILER_SHINY
	Shiny::ProfileZone pi = {NULL, Shiny::ProfileZone::STATE_HIDDEN, NULL, NULL, NULL, 0,
	                         { { 0, 0 }, { 0, 0 }, { 0, 0 } }};
	profileInformation = pi;
	profileInformation.name = pName;
	profileInformation.groups = pGroup;
	profileInformation.file = pFile;
	profileInformation.line = iLine;
	profilerCache =	&Shiny::ProfileNode::_dummy;
#endif
}

RuntimeProfileInfo::~RuntimeProfileInfo()
{
	if(bNameCopied && pName) delete[] pName;
	if(bGroupCopied && pGroup) delete[] pGroup;
	if(bFileCopied && pFile) delete[] pFile;
}

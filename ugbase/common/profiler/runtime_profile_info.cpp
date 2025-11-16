/*
 * Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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


#include "runtime_profile_info.h"
#include "string.h"

RuntimeProfileInfo::RuntimeProfileInfo(
		const char* name, bool bCopyName,
		const char* groups, bool bCopyGroup,
		const char* file, bool bCopyFile,
		int line)
: bNameCopied(false), bGroupCopied(false), bFileCopied(false),
  pName(nullptr), pGroup(nullptr), pFile(nullptr)
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
	Shiny::ProfileZone pi = {nullptr, Shiny::ProfileZone::STATE_HIDDEN, nullptr, nullptr, nullptr, 0,
	                         { { 0, 0 }, { 0, 0 }, { 0, 0 } }};
	profileInformation = pi;
	profileInformation.name = pName;
	profileInformation.groups = pGroup;
	profileInformation.file = pFile;
	profileInformation.line = iLine;
	profilerCache =	&Shiny::ProfileNode::_dummy;
#endif
#ifdef UG_PROFILER_SCOREP
	m_pHandle = SCOREP_USER_INVALID_REGION;
#endif
}

RuntimeProfileInfo::~RuntimeProfileInfo()
{
	if(bNameCopied && pName) delete[] pName;
	if(bGroupCopied && pGroup) delete[] pGroup;
	if(bFileCopied && pFile) delete[] pFile;
}

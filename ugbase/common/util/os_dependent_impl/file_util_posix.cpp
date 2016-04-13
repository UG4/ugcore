/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include <dirent.h>
#include <unistd.h>
#include <limits.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "common/util/file_util.h"
#include "common/profiler/profiler.h"
#include "common/error.h"

#ifndef PATH_MAX
// Normally, this should not happen: PATH_MAX is a POSIX constant.
// But if it happens, 1024 should be enough.
// Try to avoid this constant where possible!
#define PATH_MAX 1024
#endif

using namespace std;

namespace ug{

bool DirectoryExists(const char* dirname)
{
	DIR* curDir = opendir(dirname);
	if(!curDir)
		return false;
	
	closedir(curDir);
	return true;
}


//	This method returns a list of all directories in a directory
bool GetDirectoriesInDirectory(std::vector<std::string>& dirsOut, const char* dir)
{
	PROFILE_FUNC();
	dirsOut.clear();

	DIR* curDir = opendir(dir);

	if(!curDir) return false;

	string tFilename;

	while(dirent* entry = readdir(curDir)){
	//	get information on the file
		tFilename = dir;
		tFilename.append("/").append(entry->d_name);

		struct stat statbuf;
		stat(tFilename.c_str(), &statbuf);

		if(S_ISDIR(statbuf.st_mode))
			dirsOut.push_back(entry->d_name);
	}

	closedir(curDir);

	return true;
}

//	This method returns a list of all files in a directory
bool GetFilesInDirectory(std::vector<std::string>& filesOut, const char* dir)
{
	PROFILE_FUNC();
	filesOut.clear();

	DIR* curDir = opendir(dir);

	if(!curDir) return false;

	string tFilename;

	while(dirent* entry = readdir(curDir)){
	//	get information on the file
		tFilename = dir;
		tFilename.append("/").append(entry->d_name);

		struct stat statbuf;
		stat(tFilename.c_str(), &statbuf);

		if(S_ISREG(statbuf.st_mode))
			filesOut.push_back(entry->d_name);
	}

	closedir(curDir);

	return true;
}


bool CreateDirectoryTMP(const char *directory)
{
	return CreateDirectory(directory);
}

bool CreateDirectory(const char *directory)
{
	return mkdir(directory, 0777) == 0;
}

bool CreateDirectory(const char *directory, int mode)
{
	return mkdir(directory, mode) == 0;
}

bool CreateDirectory(std::string directory)
{
	return CreateDirectory(directory.c_str());
}

std::string GetTmpPath()
{
	return string("/tmp");
}

void ChangeDirectory(std::string dir)
{
	if(chdir(dir.c_str()) != 0){
		UG_THROW("ChangeDirectory failed.");
	}
}

std::string CurrentWorkingDirectory()
{
	char p_w_d [PATH_MAX];

	if (getcwd (p_w_d, sizeof (p_w_d)) == NULL)
		UG_THROW ("CurrentWorkingDirectory: Failed to get the current working path!");
	return std::string (p_w_d);
}

}// end of namespace

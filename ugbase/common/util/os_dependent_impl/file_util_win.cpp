/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#define _WIN32_IE 0x0501 //0x0400 //required so that SHGetSpecialFolderPath is exposed by shlobj.h

#include <windows.h>
#include <direct.h>
#include <stdlib.h>
#include <stdio.h>

#include "shlobj.h"
#include "common/util/file_util.h"
#include "common/error.h"
#include "common/log.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug {

size_t FileSize( const char *filename )
{
//	see 'https://www.securecoding.cert.org/confluence/display/c/FIO19-C.+Do+not+use+fseek%28%29+and+ftell%28%29+to+compute+the+size+of+a+regular+file'
//	for more information on the chosen implementation.
	PROFILE_FUNC();

	HANDLE file;
	LARGE_INTEGER fileSize;

	file = CreateFile(TEXT(filename), GENERIC_READ, FILE_SHARE_READ, nullptr,
					  OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, nullptr);
	UG_COND_THROW(INVALID_HANDLE_VALUE == file, "The file " << filename
	              << " could not be opened. Error " << GetLastError() );
	UG_COND_THROW(!GetFileSizeEx(file, &fileSize),
	              "Couldn't compute file size of file " << filename);
	CloseHandle(file);
	return static_cast<size_t>(fileSize.QuadPart);
}

bool DirectoryExists(const char* dirname)
{
	WIN32_FIND_DATA findData;

//	check if directory exists
	if(FindFirstFile(dirname, &findData) ==
		INVALID_HANDLE_VALUE)
	{
	//	no
		return false;
	}
	return true;
}
	
//	This method returns a list of all directories in a directory
bool GetDirectoriesInDirectory(std::vector<std::string>& dirsOut, const char* dir)
{
	dirsOut.clear();

	string expr = dir;
	expr.append("\\*");

	WIN32_FIND_DATA	findData;
	HANDLE hFind;
	
	hFind = FindFirstFile(expr.c_str(), &findData);
	if(hFind == INVALID_HANDLE_VALUE){
		return true;
	}
	
	do
	{
		if(findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY)
		{
			dirsOut.push_back(string(findData.cFileName));
		}
	}
	while(FindNextFile(hFind, &findData) != 0);
	
	FindClose(hFind);

	return true;
}

//	This method returns a list of all files in a directory
bool GetFilesInDirectory(std::vector<std::string>& filesOut, const char* dir)
{

	filesOut.clear();

	string expr = dir;
	expr.append("\\*");

	WIN32_FIND_DATA	findData;
	HANDLE hFind;
	
	hFind = FindFirstFile(expr.c_str(), &findData);
	if(hFind == INVALID_HANDLE_VALUE){
		return true;
	}
	
	do
	{
		if(!(findData.dwFileAttributes & FILE_ATTRIBUTE_DIRECTORY))
		{
			filesOut.push_back(string(findData.cFileName));
		}
	}
	while(FindNextFile(hFind, &findData) != 0);
	
	FindClose(hFind);

	return true;
}

bool CreateDirectoryTMP(const char *directory)
{
	return _mkdir(directory) == 0;
}

bool CreateDirectory(const char *directory)
{
	return _mkdir(directory) == 0;
}

bool CreateDirectory(const char *directory, int mode)
{
	return _mkdir(directory) == 0;
}

std::string GetTmpPath()
{
	static string userDataPath("");

	if(userDataPath.empty())
	{
		char tFolderPath[MAX_PATH];
		ZeroMemory(tFolderPath, MAX_PATH);

	//	init the data path
		SHGetSpecialFolderPath(0, tFolderPath, CSIDL_APPDATA, 1);

		userDataPath = tFolderPath;
	}

	return userDataPath.c_str();
}

void ChangeDirectory(std::string dir)
{
	if(chdir(dir.c_str()) != 0){
		UG_THROW("ChangeDirectory failed.");
	}
}

std::string CurrentWorkingDirectory()
{
	char * p_w_d;
	
	if ((p_w_d = _getcwd (nullptr, 0)) == nullptr )
		UG_THROW ("CurrentWorkingDirectory: Failed to get the current working path!");
	
	std::string the_pwd (p_w_d);
	free (p_w_d);
	
	return the_pwd;
}

}// end of namespace

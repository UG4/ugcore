// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)
 
#include <windows.h>

#include "common/util/file_util.h"

using namespace std;

namespace ug{

//	This method returns a list of all directories in a directory
bool GetDirectoriesInDirectory(std::vector<std::string>& dirsOut, const char* dir)
{
	dirsOut.clear();

	WIN32_FIND_DATA	findData;
	HANDLE hFind;
	
	hFind = FindFirstFile(dir, &findData);
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

	WIN32_FIND_DATA	findData;
	HANDLE hFind;
	
	hFind = FindFirstFile(dir, &findData);
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

}// end of namespace

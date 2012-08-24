// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)
 
#include <dirent.h>
#include <sys/types.h>
#include <sys/stat.h>
#include "common/util/file_util.h"
#include "common/profiler/profiler.h"

using namespace std;

namespace ug{

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


bool CreateDirectory(const char *directory, int mode)
{
	return mkdir(directory, mode) == 0;
}

std::string GetTmpPath()
{
	return string("/tmp");
}

}// end of namespace

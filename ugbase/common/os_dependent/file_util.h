// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)

#ifndef __H__UG__directory_util__
#define __H__UG__directory_util__

#include <string>
#include <vector>
#include "../ug_config.h"

namespace ug
{

///	Returns a list of all directories in a directory
UG_API bool GetDirectoriesInDirectory(std::vector<std::string>& dirsOut,
								const char* dir);

///	Returns a list of all files in a directory
UG_API bool GetFilesInDirectory(std::vector<std::string>& filesOut, const char* dir);

///	Returns true, if the specified file exists, false if not
UG_API bool FileExists(const char* filename);

}//	end of namespace

#endif

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 14.09.2011 (m,d,y)

#ifndef __H__UG__directory_util__
#define __H__UG__directory_util__

#include <string>
#include <vector>

namespace ug
{

///	This method returns a list of all directories in a directory
bool GetDirectoriesInDirectory(std::vector<std::string>& dirsOut,
								const char* dir);

///	This method returns a list of all files in a directory
bool GetFilesInDirectory(std::vector<std::string>& filesOut, const char* dir);

}//	end of namespace

#endif

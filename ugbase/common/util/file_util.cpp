/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Torbjörn Klatt
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

/**
 * \file file_util.cpp
 * \date 2012-05-15
 * \brief Implementation of OS independent file utility functions
 */

#include "common/util/file_util.h"
#include "common/ug_config.h"
#include "common/error.h"
#include "common/assert.h"
#include "common/profiler/profiler.h"
#include "path_provider.h"
#include <vector>
#include <cstring>

using namespace std;

namespace ug
{
/// !!! Serial i/o version !!!
UG_API bool FileExists( const char *filename )
{
  PROFILE_FUNC(); // since i/o

  ifstream in( filename );
  if( in.good() ) {
    in.close();
    return true;
  }
  return false;
}

/// !!! Serial i/o version !!!
// UG_API size_t FileSize( const char *filename )
// {
//   PROFILE_FUNC(); // since i/o

//   if( !FileExists( filename ) ) {
//     UG_THROW( "The file " << filename << " could not be found." );
//   }

//   ifstream ifs;
//   ifs.open( filename, ios_base::in );

//   if( !ifs.good() ) {
//     UG_THROW( "The file " << filename << " could not be opened." );
//   }
//   ifs.seekg( 0, ios_base::end );
//   size_t length = ifs.tellg();
//   ifs.close();
//   return length;
// }

bool FileTypeIs( const char* filename, const char* extension )
{
	PROFILE_FUNC();
	std::string name( filename );
	size_t iExtPos = name.find_last_of(".");
	return ( iExtPos != std::string::npos && name.substr(iExtPos).compare(extension) == 0 );
}

/// !!! Serial i/o version !!!
UG_API bool FileCompare( const char *file1, const char *file2 )
{
  PROFILE_FUNC(); // since i/o

  // Make sure, both files do exist
  if( !FileExists( file1 ) || !FileExists( file2 ) ) {
    UG_THROW( "One or both files could not be found:" << file1 << ", " << file2 );
  }

  // If both files are one and the same, we don't need to compare them any more
  if( strcmp( file1, file2 ) == 0 ) {
    return true;
  }

  // If the files have different size, we know that they are different and can stop
  size_t size1 = FileSize( file1 ), size2 = FileSize( file2 );
  if( size1 != size2 ) {
    return false;
  }
  // If both files have zero size, we can stop as well
  if( size1 == 0 && size2 == 0 ) {
    return true;
  }

  // Open the files
  ifstream f1, f2;
  f1.open( file1, ios_base::in );
  f2.open( file2, ios_base::in );

  // This should be obsolete, as it is checked by FileExists, but we want to be
  // sure.
  if( !f1.good() || !f2.good() ) {
    UG_THROW( "One or both files could not be opened:" << file1 << ", " << file2 );
  }

  string line1 = "", line2 = "";
  bool diff = false;
  while( !f1.eof() || !f2.eof() ) {
    getline( f1, line1 );
    getline( f2, line2 );
    if( ( line1 != line2 ) ||
        ( f1.eof() && !f2.eof() ) ||
        ( !f1.eof() && f2.eof() ) ) {
      diff = true;
      break;
    }
  }

  f1.close();
  f2.close();

  return !diff;
}

/// !!! Serial i/o version !!!
/// in parallel, see ParallelReadFile.
bool ReadFile(const char* filename, vector<char> &file, bool bText)
{
  PROFILE_FUNC();
  
  size_t fileSize = FileSize(filename);

  FILE *f = fopen(filename, "rb");
  if(f==NULL) return false;


  if(bText)
    file.resize(fileSize + 1);
  else
    file.resize(fileSize);

  size_t readSize = fread(&file[0], 1, fileSize, f);
  UG_COND_THROW(readSize != fileSize,
  			  "Read mismatch in ReadFile. Wrong number of bytes read: "
          << readSize << ", expected: " << fileSize);

  fclose(f);
  if(bText) file[fileSize]=0x00;
  return true;
}

/// !!! Serial i/o version !!!
string MakeTmpFile(string filename, const string &extension, bool &bSuccess)
{
	PROFILE_FUNC();  // since i/o
	bSuccess = true;
	string t = filename+extension;
	const char *name = t.c_str();
	if(!FileExists(name)) return name;
	for(int i=0;i<999999; i++)
	{
		stringstream ss;
		ss << filename << i << extension;
		name = ss.str().c_str();
		if(!FileExists(name)) return name;
	}	
	bSuccess = false;
	return "";
}

////////////////////////////////////////////////////////////////////////////////
std::string FindFileInStandardPaths(const char* filename)
{
//	first check whether the file can be loaded (e.g. absolute path or relative to woring-directory)
	std::string filenameOut = filename;
	if(FileExists(filenameOut.c_str()))
		return filenameOut;


//	Now check whether the file was specified relative to the current
//	scripting-directory
	filenameOut = PathProvider::get_current_path();
	filenameOut.append("/").append(filename);

	if(FileExists(filenameOut.c_str()))
		return filenameOut;

//	filename couldn't be located
	return "";
}


std::string FindDirInStandardPaths(const char* dirname)
{
//	first check whether the file can be loaded (e.g. absolute path or relative to woring-directory)
	std::string dirNameOut = dirname;
	if (DirectoryExists(dirNameOut.c_str()))
		return dirNameOut;

//	Now check whether the file was specified relative to the current
//	scripting-directory
	dirNameOut = PathProvider::get_current_path();
	dirNameOut.append("/").append(dirname);

	if(DirectoryExists(dirNameOut.c_str()))
		return dirNameOut;

//	filename couldn't be located
	return "";
}



} // namespace ug

// EOF

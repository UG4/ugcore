/*
 * Copyright (c) 2011-2014:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter, Torbjörn Klatt
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
 * \file file_util.h
 * \date 2012-05-15
 * \brief File utility functions
 * \details Utility function, which depend on functionality specific to
 * operating systems are implemented in common/os_dependent/file_util_*.cpp
 */

#ifndef __H__UG__FILE_UTIL__
#define __H__UG__FILE_UTIL__

#include <fstream>
#include <vector>
#include <string>
#include "common/ug_config.h"

namespace ug
{

/// \addtogroup ugbase_common_io
/// \{

/**
 * \brief Returns a list of all directories in a directory
 *
 * The returned list contains as well the directory scanned as <tt>.</tt> and its
 * parent directory as <tt>..</tt>.
 *
 * \param[out]  dirsOut string vector for holding the list of directories
 * \param[in]   dir     path and name of the directory to scan
 *
 * \note The implementation relies on OS specific instructions.
 *       Implementations for POSIX-UNIX and Windows are available.
 */
UG_API bool GetDirectoriesInDirectory(std::vector<std::string>& dirsOut, const char* dir);

/**
 * \brief Returns a list of all files in a directory
 *
 * \param[out]  filesOut  string vector for holding the list of files
 * \param[in]   dir       path and name of the directory to scan
 *
 * \note The implementation relies on OS specific instructions.
 *       Implementations for POSIX-UNIX and Windows are available.
 */
UG_API bool GetFilesInDirectory(std::vector<std::string>& filesOut, const char* dir);

/**
 * \brief Checks the existence of a given directory
 *
 * \param[in] name of the directory to be checked
 * \return \c true if the specified directory exists, \c false otherwise
 */
UG_API bool DirectoryExists(const char* dirname);
static inline bool DirectoryExists( std::string filename)
{
	return DirectoryExists(filename.c_str());
}
/**
 * \brief Checks the existence of a given file
 *
 * \param[in] filename path and name of the file to be checked
 * \return \c true if the specified file exists, \c false otherwise
 */
UG_API bool FileExists( const char *filename );
static inline bool FileExists( std::string filename)
{
	return FileExists(filename.c_str());
}

/**
 * \brief Returns size of the specified file in bytes
 *
 * \param[in] filename Path and name of the file
 * \return Size of \c filename in bytes
 * \throws std::runtime_error if the specified file does not exist or could not
 *                            be opened for reading.
 */
UG_API size_t FileSize( const char *filename );
static inline size_t FileSize( std::string filename)
{
	return FileSize(filename.c_str());
}

/**
 * \brief Check filename extension
 * 
 * \param[in] filename path and anme of the file
 * \param[in] extension desired file extension (including leading dot)
 * \return \c true if \c filename has given extension, \c false otherwise
 */
UG_API bool FileTypeIs( const char *filename, const char *extension );

/**
 * \brief Creates a directory
 *
 * \param[in] directory name of the directory
 * \param[in] mode	(optional, default 0777) Sets ownership options for the file.
 *					Ignored on windows.
 * \return true if successfull
 * \{
 */
UG_API bool CreateDirectoryTMP(const char *directory);
UG_API bool CreateDirectory(const char *directory);
UG_API bool CreateDirectory(const char *directory, int mode);
static inline bool CreateDirectory(std::string directory)
{
	return CreateDirectoryTMP(directory.c_str());
}
/** \} */


/**
 * \brief Compares two files by their content
 *
 * Comparison is first done on file size using ug::FileSize.
 * If the two files have same size, then the comparison is done line wise based
 * on std::string comparison.
 *
 * \param[in] file1 Path and name of the first file
 * \param[in] file2 Path and name of the second file
 * \return \c true if they are the same one or different ones but with same
 *         content; \c false otherwise
 * \throws std::runtim_error if one or both files can not be found or opened for
 *                           reading.
 */
UG_API bool FileCompare( const char *file1, const char *file2 );

/**
 * \brief Returns a path to which an application may write temporary files.
 *
 * On unix systems, this is normally "/tmp". On windows the users AppData path
 * or something similar is returned.
 */
UG_API std::string GetTmpPath();

/**
 * \param filename filename to read
 * \param file vector to put whole file into
 * \param bText if true, open file with r, otherwise rb
 * used in \sa ParallelReadFile and Script
 */
UG_API bool ReadFile(const char* filename, std::vector<char> &file, bool bText);


/**
 * \param filename filename including path (e.g. /Users/horst/file1)
 * \param extension the extension to be added to the file e.g. txt
 * \param bSuccess true if successful, false otherwise
 * \return filenameXXXXXX.extension, where XXXXXX is some number between 1 and 999999,
 *			so that the file doesn't exist yet.
 */
UG_API std::string MakeTmpFile(std::string filename, const std::string &extension,
		bool &bSuccess);

///	Changes the current directory
UG_API void ChangeDirectory(std::string dir);


////////////////////////////////////////////////////////////////////////////////
///	searches the file in the standard paths.
/**	The checking order is the following:
 *		- absolute or relative path to working directory
 * 		- relative path regarding the current script path (PathProvider::get_current_path())
 *
 * \return	the full path (including the filename) at which the file was found, or an
 *			empty string if no matching file was found.
 */
std::string FindFileInStandardPaths(const char* filename);

////////////////////////////////////////////////////////////////////////////////
///	searches the directory in the standard paths.
/**	The checking order is the following:
 *		- absolute or relative path to working directory
 * 		- relative path regarding the current script path (PathProvider::get_current_path())
 *
 * \return	the full path (including the dir name) at which the dir was found, or an
 *			empty string if no matching file was found.
 */
std::string FindDirInStandardPaths(const char* dirname);

/** Current working directory
 *
 * \return	the current working directory as a string
 */
std::string CurrentWorkingDirectory();

// end group ugbase_common_io
/// \}

} // namespace ug

#endif // __H__UG__FILE_UTIL__

// EOF

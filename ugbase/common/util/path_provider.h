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

#ifndef __H__UG__path_provider__
#define __H__UG__path_provider__

#include <string>
#include <map>
#include <stack>

#include "common/util/file_util.h"
#include "common/util/os_info.h"


namespace ug {

/// \addtogroup ugbase_common_util
/// \{

////////////////////////////////////////////////////////////////////////
///	Constants used by PathProvider
enum PathTypes
{
	BIN_PATH = 0,	///< path in which the binary lies
	SCRIPT_PATH,
	ROOT_PATH,
	PLUGIN_PATH,
	APPS_PATH,		///< path in which the application-scripts lie

//	always last
	MAX_PATH_CONSTANT
};

////////////////////////////////////////////////////////////////////////
///	Singleton which stores common paths and a stack of current paths.
/**	All paths are initially set to "".
 *
 * Note that all public methods of PathProvider are static. That means
 * you have to call them through the :: operator. E.g.
 *
 * \code
 * std::string appPath = PathProvider::get_path(BIN_PATH);
 * \endcode
 */
class PathProvider
{
	public:
	///	sets the path for the given constant.
	/**
	 * \param pathType	should be one of the constants enumerated in PathTypes
	 * 					or a used defined constant starting from
	 * 					MAX_PATH_CONSTANT + 1.
	 */
		static inline void set_path(PathTypes pathType, const std::string& path)
		{inst().m_map[pathType] = path;}

	///	returns the path associated with the given constant.
	/**
	 * \param pathType	should be one of the constants enumerated in PathTypes
	 * 					or a used defined constant starting from
	 * 					MAX_PATH_CONSTANT + 1.
	 */
		static inline const std::string& get_path(PathTypes pathType)
		{return inst().m_map[pathType];}

	///	returns true, if the path associated with the given constant exists.
	/**
	 * \param pathType	should be one of the constants enumerated in PathTypes
	 * 					or a used defined constant starting from
	 * 					MAX_PATH_CONSTANT + 1.
	 */
		static inline bool has_path(PathTypes pathType)
		{return inst().m_map.find(pathType) != inst().m_map.end();}

	///	returns the current path
	/**	current paths are stored in a stack. The top of the stack is considered
	 * to be the most current path and is returned by this method.
	 *
	 * \param defPath	(optional) If the stack is empty, the path associated with
	 * 					defPath is returned. By default defPath is set to BIN_PATH.
	 */
		static inline const std::string& get_current_path(PathTypes defPath = BIN_PATH)
		{
			if(inst().m_curPaths.empty())
				return get_path(defPath);
			return inst().m_curPaths.top();
		}

	///	returns true if a current path exists, false if not.
		static inline bool has_current_path()
		{return !inst().m_curPaths.empty();}

	///	pushes a path to the stack of current paths
		static inline void push_current_path(const std::string& path)
		{inst().m_curPaths.push(path);}

	///	pops a path from the stack of current paths
		static inline void pop_current_path()
		{inst().m_curPaths.pop();}

	///	clears the stack of current paths. This makes sense if an error was catched.
		static inline void clear_current_path_stack()
		{while(has_current_path()) {pop_current_path();}}

	/**
	 * @param relativeFilename (in) relative filename
	 * @param absoluteFilename (out) absolute filename
	 * @return true if file exists relative to current path
	 */
		static inline bool get_filename_relative_to_current_path(const std::string &relativeFilename, std::string &absoluteFilename)
		{
			const char* pathSep = GetPathSeparator();

			if(has_current_path() == false) return false;
			absoluteFilename = get_current_path() + pathSep + relativeFilename;
			return FileExists(absoluteFilename.c_str());
		}


	/**
	 * @param relativeDirname (in) relative directory name
	 * @param absoluteDirname (out) absolute directory name
	 * @return true if directory exists relative to current path
	 */
		static inline bool get_dirname_relative_to_current_path(const std::string &relativeDirname, std::string &absoluteDirname)
		{
			const char* pathSep = GetPathSeparator();

			if (!has_current_path())
				return false;
			absoluteDirname = get_current_path() + pathSep + relativeDirname;
			return DirectoryExists(absoluteDirname.c_str());
		}

	/**
	 * @param relativeFilename (in) relative filename
	 * @param absoluteFilename (out) absolute filename
	 * @return true if file exists relative to current path
	 */
		static inline bool get_filename_relative_to_path(PathTypes pathType, const std::string &relativeFilename, std::string &absoluteFilename)
		{
			const char* pathSep = GetPathSeparator();

			if(has_path(pathType) == false) return false;
			absoluteFilename = get_path(pathType) + pathSep + relativeFilename;
			return FileExists(absoluteFilename.c_str());
		}

	/**
	 * @param relativeDirname (in) relative directory name
	 * @param absoluteDirname (out) absolute directory name
	 * @return true if directory exists relative to current path
	 */
		static inline bool get_dirname_relative_to_path(PathTypes pathType, const std::string &relativeDirname, std::string &absoluteDirname)
		{
			const char* pathSep = GetPathSeparator();

			if (!has_path(pathType))
				return false;
			absoluteDirname = get_path(pathType) + pathSep + relativeDirname;
			return DirectoryExists(absoluteDirname.c_str());
		}
	private:
		PathProvider() = default;
		PathProvider(const PathProvider&)	{}

		static PathProvider& inst()
		{
			static PathProvider pp;
			return pp;
		}

	private:
		std::map<PathTypes, std::string>	m_map;
		std::stack<std::string>		m_curPaths;
};

// end group ugbase_common_util
/// \}

}//	end of namespace

#endif

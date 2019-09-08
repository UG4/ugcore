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

#include "../util_overloaded.h"
#include "ug.h"
#include "bridge/bridge.h"
#include "common/stopwatch.h"
#include "common/util/file_util.h"
#include "common/util/path_provider.h"
#include "common/util/table.h"
#include "common/util/variant.h"
#include "registry/registry.h"
#if defined (__APPLE__) || defined (__linux__)
	#include "common/util/mem_info.h"
#endif

#include <cstdlib>
#include <string>

using namespace std;

namespace ug{

// determines when to show the progress bar
extern size_t g_minSecondUntilProgress;

namespace bridge{

/// \defgroup util_bridge Utility Bridge
/// \ingroup misc_bridge
/// \{

static string GetRootPath()
{return PathProvider::get_path(ROOT_PATH);}

static string GetAppsPath()
{return PathProvider::get_path(APPS_PATH);}

static string GetBinPath()
{return PathProvider::get_path(BIN_PATH);}

static string GetScriptPath()
{return PathProvider::get_path(SCRIPT_PATH);}

static string GetCurrentPath()
{return PathProvider::get_current_path();}

/**
 * This function executes the given string in a system shell (if available)
 *
 * \param[in]	cmd		The shell command to execute
 * \returns 	1 		if no shell available
 */
static int ExecuteSystemCommand(const char* cmd)
{
//	check that some shell exists
	if(!system(NULL)) return 1;

//	run the command
	return system(cmd);
}

void int_srand(int seed)
{
	srand((unsigned int)seed);
}

void SetMinSecondsUntilProgress(size_t s)
{
	g_minSecondUntilProgress=s;
}

std::vector<std::string> GetFilesInDir(const char* dir)
{
	std::vector<std::string> files;
	GetFilesInDirectory(files, dir);
	return files;
}

std::vector<std::string> GetDirsInDir(const char* dir)
{
	std::vector<std::string> dirs;
	GetDirectoriesInDirectory(dirs, dir);
	return dirs;
}

void RegisterBridge_Util(Registry& reg, string parentGroup)
{
	string grp(parentGroup);
	grp.append("/Util");

	reg.add_function("ug_get_root_path", &GetRootPath, grp,
	                 "pathName", "", "Returns ug's root path");

	reg.add_function("ug_get_apps_path", &GetAppsPath, grp,
	                 "pathName", "", "Returns the path in which ug's apps are stored");

	reg.add_function("ug_get_bin_path", &GetBinPath, grp,
	                 "pathName", "", "Returns the path in which the ug executable lies");

	reg.add_function("ug_get_script_path", &GetScriptPath, grp,
	                 "pathName", "", "Returns the script path");

	reg.add_function("ug_get_current_path", &GetCurrentPath, grp,
	                 "pathName", "", "Returns the current path");

	reg.add_function("ug_set_root_path", &SetRootPath, grp,
	                 "", "pathName", "Sets the paths relative to passed root path");

	reg.add_function("ExecuteSystemCommand", &ExecuteSystemCommand, grp,
	                 "success", "command", "Executes a command in the system shell");

	reg.add_function("srand", int_srand, grp, 
	                 "", "seed", "The pseudo-random number generator is initialized using the argument passed as seed.")
		.add_function("ug_file_exists", OVERLOADED_FUNCTION_PTR(bool, FileExists, (const char*)), grp,
	                  "exists", "", "Returns true if a path exists, false if not.")
		.add_function("exit", &UGForceExit, grp,
	                  "", "", "Immediatly terminates the application.")
		.add_function("quit", &UGForceExit, grp,
	                  "", "", "Immediatly terminates the application.");

	reg.add_function("CreateDirectory", static_cast<bool (*)(const char *) > (&CreateDirectoryTMP) );
	reg.add_function("DirectoryExists", static_cast<bool (*)(const char *) > (&DirectoryExists) );
	reg.add_function("FileExists", static_cast<bool (*)(const char *) > (&FileExists) );
	reg.add_function("GetFilesInDir", GetFilesInDir);
	reg.add_function("GetDirsInDir", GetDirsInDir);
	reg.add_function("GetTmpPath", GetTmpPath);
	reg.add_function("ChangeDirectory", ChangeDirectory);
	reg.add_function("CurrentWorkingDirectory", CurrentWorkingDirectory);
	reg.add_function("FileCompare", FileCompare);
	reg.add_function("SetMinSecondsUntilProgress", SetMinSecondsUntilProgress, grp, "", "seconds", "determines after which time a progress bar can show up");

	reg.add_function("FindFileInStandardPaths", FindFileInStandardPaths);

	{
		typedef Variant T;
		reg.add_class_<T>("Variant", grp)
		.add_constructor()
		.add_constructor<void (*)(bool)>()
		.add_constructor<void (*)(int)>()
		.add_constructor<void (*)(size_t)>()
		.add_constructor<void (*)(float)>()
		.add_constructor<void (*)(double)>()
		.add_constructor<void (*)(const char*)>()
		.add_constructor<void (*)(const Variant&)>()
		.add_method("to_bool", &T::to_bool)
		.add_method("to_int", &T::to_int)
		.add_method("to_size_t", &T::to_size_t)
		.add_method("to_float", &T::to_float)
		.add_method("to_double", &T::to_double)
		.add_method("to_number", &T::to_number)
		.add_method("to_string", &T::to_c_string)
		.construct_as_smart_pointer();
	}

	{
		typedef StringTable T;
		reg.add_class_<T>("StringTable", grp)
		.add_constructor()
		.add_constructor< void (*) ( size_t numRows, size_t numCols) > ()
		.add_method("set", &T::set)
		.add_method("get", &T::get)
		.add_method("add_rows", &T::add_rows)
		.add_method("add_cols", &T::add_cols)
		.add_method("num_rows", &T::num_rows)
		.add_method("num_cols", &T::num_cols)
		.add_method("set_col_alignment", &T::set_col_alignment)
		.add_method("set_col_alignments", &T::set_col_alignments)
		.add_method("set_default_col_alignment", &T::set_default_col_alignment)

		.add_method("set_row_seperator", &T::set_row_seperator)
		.add_method("set_row_seperators", &T::set_row_seperators)
		.add_method("set_default_row_seperator", &T::set_default_row_seperator)

		.add_method("set_col_seperator", &T::set_col_seperator)
		.add_method("set_col_seperators", &T::set_col_seperators)
		.add_method("set_default_col_seperator", &T::set_default_col_seperator)

		.add_method("to_latex", &T::to_latex)
		.add_method("to_csv", &T::to_csv)
		.add_method("__tostring", &T::to_string)

		.add_method("transpose", &T::transpose)
		
		.construct_as_smart_pointer();
	}

	{
		// Matlab like stop watch
		 typedef CuckooClock T;
		 reg.add_class_<T>("CuckooClock", grp)
			.add_constructor()
		    .add_method("tic", &T::tic)
		    .add_method("toc", &T::toc)
		    .add_method("cuckoo", &T::cuckoo);
	}

#if defined (__APPLE__) || defined (__linux__)
	// MemInfo provides information about memory usage
	{
		typedef MemInfo T;
		string name = string("MemInfo");
		reg.add_class_<T>(name, grp)
			.add_constructor()
			.add_method("memory_consumption", &T::memory_consumption, "", "", "")
			.add_method("local_resident_memory", &T::local_resident_memory, "", "", "")
			.add_method("local_virtual_memory", &T::local_virtual_memory, "", "", "")
			.add_method("global_resident_memory", &T::global_resident_memory, "", "", "")
			.add_method("global_virtual_memory", &T::global_virtual_memory, "", "", "")
			.add_method("max_resident_memory", &T::max_resident_memory, "", "", "")
			.add_method("max_virtual_memory", &T::max_virtual_memory, "", "", "")

			.set_construct_as_smart_pointer(true);
	}
#endif

}

// end group util_bridge
/// \}

}// end of namespace bridge
}// end of namespace ug

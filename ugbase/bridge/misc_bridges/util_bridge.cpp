// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)

#include "common/util/path_provider.h"
#include "registry/registry.h"
#include "bridge/bridge.h"
#include "common/util/file_util.h"
#include "../util_overloaded.h"
#include "ug.h"
#include "common/util/table.h"
#include "common/stopwatch.h"

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

static string GetGridPath()
{return PathProvider::get_path(GRID_PATH);}

static string GetAppsPath()
{return PathProvider::get_path(APPS_PATH);}

static string GetBinPath()
{return PathProvider::get_path(BIN_PATH);}

static string GetDataPath()
{return PathProvider::get_path(DATA_PATH);}

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

	reg.add_function("ug_get_grid_path", &GetGridPath, grp,
	                 "pathName", "", "Returns the path in which ug's standard grids are stored.");

	reg.add_function("ug_get_apps_path", &GetAppsPath, grp,
	                 "pathName", "", "Returns the path in which ug's apps are stored");

	reg.add_function("ug_get_bin_path", &GetBinPath, grp,
	                 "pathName", "", "Returns the path in which the ug executable lies");

	reg.add_function("ug_get_data_path", &GetDataPath, grp,
	                 "pathName", "", "Returns the data path");

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
	reg.add_function("FileCompare", FileCompare);
	reg.add_function("SetMinSecondsUntilProgress", SetMinSecondsUntilProgress, grp, "", "seconds", "determines after which time a progress bar can show up");

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

}

// end group util_bridge
/// \}

}// end of namespace bridge
}// end of namespace ug

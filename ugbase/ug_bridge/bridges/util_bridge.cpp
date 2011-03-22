// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 22.03.2011 (m,d,y)

#include "common/util/path_provider.h"
#include "../registry.h"
#include "../ug_bridge.h"

namespace ug{
namespace bridge{

static std::string GetAppPath()
{return PathProvider::get_path(APP_PATH);}

static std::string GetDataPath()
{return PathProvider::get_path(DATA_PATH);}

static std::string GetScriptPath()
{return PathProvider::get_path(SCRIPT_PATH);}

static std::string GetCurrentPath()
{return PathProvider::get_current_path();}



bool RegisterUtilInterface(Registry& reg, const char* parentGroup)
{
	std::string grpStr(parentGroup);
	grpStr.append("/util");

	reg.add_function("ug_get_app_path", &GetAppPath, grpStr.c_str(),
					 "pathName", "", "Returns the application path");

	reg.add_function("ug_get_data_path", &GetDataPath, grpStr.c_str(),
					 "pathName", "", "Returns the data path");

	reg.add_function("ug_get_script_path", &GetScriptPath, grpStr.c_str(),
					 "pathName", "", "Returns the script path");

	reg.add_function("ug_get_current_path", &GetCurrentPath, grpStr.c_str(),
					 "pathName", "", "Returns the current path");

	return true;
}

}// end of namespace bridge
}// end of namespace ug

// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Nov 25, 2013

#include <string>
#include <cstring>
#include "common/log.h"
#include "bridge/bridge.h"
#include "common/profiler/profiler.h"
#include "path_provider.h"
#include "plugin_util.h"
#include "embedded_plugins.h"	// created by cmake in the build directory

using namespace std;

namespace ug{

bool PluginLoaded(const string &name)
{
	string searchStr;
	searchStr.append(" ").append(name).append(" ");
	return (strstr(ListOfEmbeddedPlugins(), searchStr.c_str()) != NULL);
}

bool LoadPlugins(const char*, string parentGroup, bridge::Registry& reg)
{
	InitializeEmbeddedPlugins(&reg, parentGroup);
	return true;
}

bool UnloadPlugins()
{
	return true;
}

}// end of namespace

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#ifndef __H__UG__plugin_util__
#define __H__UG__plugin_util__

#include <string>
#include "common/ug_config.h"

namespace ug
{

///	Loads all plugins in the given path.
UG_API bool LoadPlugins(const char* pluginPath, std::string parentGroup);
UG_API bool UnloadPlugins();
UG_API bool PluginLoaded(const std::string &name);

}//	end of namespace

#endif

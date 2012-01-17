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
bool UG_API LoadPlugins(const char* pluginPath, std::string parentGroup);

}//	end of namespace

#endif

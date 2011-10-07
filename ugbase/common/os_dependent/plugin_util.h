// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#ifndef __H__UG__plugin_util__
#define __H__UG__plugin_util__

#include <string>

namespace ug
{

///	Loads all plugins in the given path.
bool LoadPlugins(const char* pluginPath, std::string);

}//	end of namespace

#endif

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#ifndef __H__UG__plugin_util__
#define __H__UG__plugin_util__

namespace ug
{

///	Loads the specified dynamic library and runs init_ug_plugin() on it.
bool LoadPlugin(const char* pluginName);

}//	end of namespace

#endif

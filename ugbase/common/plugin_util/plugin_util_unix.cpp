// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)
 
#include <dlfcn.h>
#include "common/log.h"

namespace ug{

bool LoadPlugin(const char* pluginName)
{
//	Thanks to Michael Hoffer for his examples!

//	load the library
	void* libHandle = dlopen(pluginName, RTLD_LAZY);
	if(!libHandle){
		UG_LOG("Failed to load plugin " << pluginName << ". File not found.\n");
		return false;
	}

//	find the init_ug_plugin function
	typedef void (*FctInitPlugin)();

	FctInitPlugin fctInitPlugin = (FctInitPlugin) dlsym(libHandle, "init_ug_plugin");

	if(!fctInitPlugin){
		UG_LOG("Method 'void init_ug_plugin()' could not be found in " << pluginName << ".\n");
		return false;
	}

//	call the init method
	fctInitPlugin();

//	we're done in here.
	return true;
}

}// end of namespace

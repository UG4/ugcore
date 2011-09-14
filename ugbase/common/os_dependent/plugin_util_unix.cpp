// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)
 
#include <dlfcn.h>
#include <string>
#include "common/log.h"
#include "common/util/path_provider.h"
#include "ug_bridge/ug_bridge.h"
#include "common/os_dependent/file_util.h"

using namespace std;

namespace ug{

bool LoadPlugins(const char* pluginPath)
{
	typedef void (*FctInitPlugin)(ug::bridge::Registry*);

//	first we'll try to find all plugins in the given path
	vector<string> files;

	GetFilesInDirectory(files, pluginPath);

	//UG_LOG("Loading plugins (from " << pluginPath << "):");

	bridge::Registry& reg = bridge::GetUGRegistry();

	for(size_t i = 0; i < files.size(); ++i){
		//UG_LOG(" " << files[i]);

		string fullPluginName(pluginPath);
		fullPluginName.append("/").append(files[i]);

		void* libHandle = dlopen(fullPluginName.c_str(), RTLD_LAZY);

		if(!libHandle){
			UG_LOG("(failed)");
			continue;
		}

		std::string fctName("InitUGPlugin");
		//fctName.append("_").append(pluginName); //pluginName not yet known...

	//	find the init_ug_plugin function
		FctInitPlugin fctInitPlugin = (FctInitPlugin) dlsym(libHandle, fctName.c_str());

		if(!fctInitPlugin){
			UG_LOG("(failed)");
			continue;
		}

	//	call the init method
		fctInitPlugin(&reg);
	}

//	make sure that the registry is updated
	reg.registry_changed();

	//UG_LOG(endl);

	return true;
}

}// end of namespace

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include <dlfcn.h>
#include <string>
#include "common/log.h"
#include "common/util/path_provider.h"
#include "bridge/bridge.h"
#include "common/os_dependent/file_util.h"

using namespace std;

namespace ug{

bool LoadPlugins(const char* pluginPath, std::string parentGroup)
{
	PROFILE_FUNC();
	typedef void (*FctInitPlugin)(ug::bridge::Registry*, std::string);

//	first we'll try to find all plugins in the given path
	vector<string> files;

	GetFilesInDirectory(files, pluginPath);

	bridge::Registry& reg = bridge::GetUGRegistry();

	for(size_t i = 0; i < files.size(); ++i){
		string fullPluginName(pluginPath);
		fullPluginName.append("/").append(files[i]);

		void* libHandle = dlopen(fullPluginName.c_str(), RTLD_LAZY);

		if(!libHandle){
			UG_LOG("PLUGIN-ERROR: Couldn't open plugin " << files[i] << endl);
			UG_LOG("Error message: " << dlerror() << endl);
			UG_LOG("NOTE: This could be due to incompatible build settings in ugshell and the plugin.\n");
			continue;
		}

		std::string fctName("InitUGPlugin");

	//	find the init_ug_plugin function
		FctInitPlugin fctInitPlugin = (FctInitPlugin) dlsym(libHandle, fctName.c_str());

		if(!fctInitPlugin){
			UG_LOG("PLUGIN-ERROR: Couldn't find entry point 'InitUGPlugin' in plugin " << files[i] << endl);
			continue;
		}

	//	call the init method
		fctInitPlugin(&reg, parentGroup);
	}

//	make sure that the registry is updated
	reg.registry_changed();

	return true;
}

}// end of namespace

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 15.09.2011 (m,d,y)

#include <windows.h>
#include <string>
#include "plugin_util.h"
#include "common/log.h"
#include "common/util/path_provider.h"
#include "bridge/bridge.h"
#include "../file_util.h"
#include "../os_info.h"

using namespace std;

namespace ug{

bool LoadPlugins(const char* pluginPath, std::string parentGroup)
{
	typedef void (*FctInitPlugin)(ug::bridge::Registry*, std::string);

	bool bSuccess = true;

//	first we'll try to find all plugins in the given path
	vector<string> files;

	GetFilesInDirectory(files, pluginPath);

	//UG_ERR_LOG("Loading plugins (from " << pluginPath << "):");

	bridge::Registry& reg = bridge::GetUGRegistry();

	for(size_t i = 0; i < files.size(); ++i){
		//UG_ERR_LOG(" " << files[i]);

		string fullPluginName(pluginPath);
		fullPluginName.append(GetPathSeparator()).append(files[i]);

		HMODULE libHandle = LoadLibrary(fullPluginName.c_str());

		if(!libHandle){
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't open plugin " << files[i]);
			UG_ERR_LOG(" (" << GetLastError() << ")\n");
			UG_ERR_LOG("NOTE: This could be due to incompatible build settings in ugshell and the plugin.\n");
			bSuccess = false;
			continue;
		}

		std::string fctName("InitUGPlugin");
		//fctName.append("_").append(pluginName); //pluginName not yet known...

	//	find the init_ug_plugin function
		FctInitPlugin fctInitPlugin =
			(FctInitPlugin) GetProcAddress(libHandle, fctName.c_str());

		if(!fctInitPlugin){
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't find entry point 'InitUGPlugin' in plugin " << files[i]);
			UG_ERR_LOG(" (" << GetLastError() << ")\n");
			bSuccess = false;
			continue;
		}

	//	call the init method
		fctInitPlugin(&reg, parentGroup);
	}

//	make sure that the registry is updated
	reg.registry_changed();

	//UG_ERR_LOG(endl);

	return bSuccess;
}


bool UnloadPlugins()
{
	// please implement
	return false;
}

}// end of namespace

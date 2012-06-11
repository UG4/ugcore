// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include <dlfcn.h>
#include <string>
#include "common/log.h"
#include "common/util/path_provider.h"
#include "bridge/bridge.h"
#include "../file_util.h"
#include "common/profiler/profiler.h"
#include "../os_info.h"

using namespace std;

namespace ug{

static vector<void*> loadedPlugins;

bool LoadPlugins(const char* pluginPath, std::string parentGroup)
{
	PROFILE_FUNC();
	typedef void (*FctInitPlugin)(ug::bridge::Registry*, std::string);

	bool bSuccess = true;

//	first we'll try to find all plugins in the given path
	vector<string> files;

	GetFilesInDirectory(files, pluginPath);

	bridge::Registry& reg = bridge::GetUGRegistry();

	int prefixLen = strlen(GetDynamicLibraryPrefix());
	int suffixLen = strlen(GetDynamicLibrarySuffix()) + 1; // include '.'

	for(size_t i = 0; i < files.size(); ++i){
	//	extract the plugins name from the file-name
		int nameStart = prefixLen;
		int nameLength = (int)files[i].size() - suffixLen - nameStart;

		if(nameLength <= 0)
		{
			UG_ERR_LOG("\nCouldn't extract plugin-name from filename " << files[i] <<
					". Ignoring plugin.\n");
			bSuccess = false;
			continue;
		}

		string pluginName = files[i].substr(nameStart, nameLength);

	//	load the library
		string fullPluginName(pluginPath);
		fullPluginName.append("/").append(files[i]);

		void* libHandle = dlopen(fullPluginName.c_str(), RTLD_LAZY);

		if(!libHandle){
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't open plugin " << files[i] << endl);
			UG_ERR_LOG("Error message: " << dlerror() << endl);
			UG_ERR_LOG("NOTE: This could be due to incompatible build settings in ugshell and the plugin.\n");
			bSuccess = false;
			continue;
		}


	//	find the plugins init function
		string fctName("InitUGPlugin_");
		fctName.append(pluginName);

		FctInitPlugin fctInitPlugin = (FctInitPlugin) dlsym(libHandle, fctName.c_str());

		if(!fctInitPlugin){
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't find entry point " <<  fctName
					<< " in plugin " << files[i] << endl);
			dlclose(libHandle);
			bSuccess = false;
			continue;
		}

	//	call the init method
		fctInitPlugin(&reg, parentGroup);
		loadedPlugins.push_back(libHandle);
	}

//	make sure that the registry is updated
	reg.registry_changed();

	return bSuccess;
}

bool UnloadPlugins()
{
	typedef void (*FctFinalizePlugin)();
	for(size_t i=0; i<loadedPlugins.size(); ++i)
	{
		std::string fctName("FinalizeUGPlugin");
		void *libHandle = loadedPlugins[i];

		FctFinalizePlugin fctFinalizePlugin = (FctFinalizePlugin) dlsym(libHandle, fctName.c_str());

		if(fctFinalizePlugin) fctFinalizePlugin();

		dlclose(libHandle);
	}
	loadedPlugins.clear();
	return true;
}

}// end of namespace

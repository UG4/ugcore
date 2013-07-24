// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include <string>
#include "common/log.h"
#include "bridge/bridge.h"
#include "common/profiler/profiler.h"
#include "path_provider.h"
#include "file_util.h"
#include "os_info.h"
#include "dynamic_library_util.h"

#include "plugin_util.h"

using namespace std;

namespace ug{

namespace {
	vector<DynLibHandle> loadedPlugins;
	vector<string> loadedPluginNames;
}

bool PluginLoaded(const string &name)
{
	return find(loadedPluginNames.begin(), loadedPluginNames.end(), name) != loadedPluginNames.end();
}

bool LoadPlugins(const char* pluginPath, string parentGroup, bridge::Registry& reg)
{
	PROFILE_FUNC();
	typedef void (*FctInitPlugin)(ug::bridge::Registry*, string);

	bool bSuccess = true;

//	first we'll try to find all plugins in the given path
	vector<string> files;

	GetFilesInDirectory(files, pluginPath);

	string prefixStr = GetDynamicLibraryPrefix();
	string suffixStr = string(".").append(GetDynamicLibrarySuffix());

	int prefixLen = prefixStr.size();
	int suffixLen = suffixStr.size(); // includes '.'

	for(size_t i = 0; i < files.size(); ++i){

	//	extract the plugins name from the file-name
		int nameStart = prefixLen;
		int nameLength = (int)files[i].size() - suffixLen - nameStart;

	//	exclude MAC OS X hidden folder custom file ".DS_Store" from plugin consideration
		#ifdef __APPLE__
			if(files[i].compare(".DS_Store") == 0) continue;
		#endif

	//	check that plugin name can exist
		if(nameLength <= 0)
		{
			UG_ERR_LOG("Plugin-filename '" << files[i] <<
					"' too short. Ignoring plugin.\n");
			bSuccess = false;
			continue;
		}

	//	check for prefix
		if(files[i].compare(0,prefixLen,prefixStr) != 0)
		{
			UG_ERR_LOG("Plugin-filename '" << files[i] << "' does not "
					"start with Plugin prefix '"<<prefixStr<<"'. Ignoring plugin.\n");
			bSuccess = false;
			continue;
		}

	//	check for suffix
		if(files[i].compare(files[i].size()-suffixLen,suffixLen,suffixStr) != 0)
		{
			UG_ERR_LOG("Plugin-filename '" << files[i] << "' does not "
					"end with Plugin suffix '"<<suffixStr<<"'. Ignoring plugin.\n");
			bSuccess = false;
			continue;
		}

	//	load the library
		string fullPluginName(pluginPath);
		fullPluginName.append("/").append(files[i]);

		DynLibHandle libHandle = OpenLibrary(fullPluginName.c_str());

		if(!libHandle){
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't open plugin " << files[i] << endl);
			UG_ERR_LOG("NOTE: This could be due to incompatible build settings in ugshell and the plugin.\n");
			bSuccess = false;
			continue;
		}


	//	find the plugins init function
		string pluginName = files[i].substr(nameStart, nameLength);
		string fctName("InitUGPlugin_");
		fctName.append(pluginName);

		FctInitPlugin fctInitPlugin =
				(FctInitPlugin) GetLibraryProcedure(libHandle, fctName.c_str());

		if(!fctInitPlugin){
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't find entry point " <<  fctName
					<< " in plugin " << files[i] << endl);
			CloseLibrary(libHandle);
			bSuccess = false;
			continue;
		}

	//	call the init method
		fctInitPlugin(&reg, parentGroup);
		loadedPlugins.push_back(libHandle);
		loadedPluginNames.push_back(pluginName);
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
		string fctName("FinalizeUGPlugin");
		DynLibHandle libHandle = loadedPlugins[i];

		FctFinalizePlugin fctFinalizePlugin =
				(FctFinalizePlugin) GetLibraryProcedure(libHandle, fctName.c_str());

		if(fctFinalizePlugin) fctFinalizePlugin();

		CloseLibrary(libHandle);
	}
	loadedPlugins.clear();
	return true;
}

}// end of namespace

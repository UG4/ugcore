// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 15.09.2011 (m,d,y)

#include <windows.h>
#include <string>
#include <vector>
#include "../plugin_util.h"
#include "common/log.h"
#include "common/util/path_provider.h"
#include "bridge/bridge.h"
#include "../file_util.h"
#include "../os_info.h"

using namespace std;

namespace ug{

static std::vector<std::string> loadedPluginNames;
bool PluginLoaded(const std::string &name)
{
	return std::find(loadedPluginNames.begin(), loadedPluginNames.end(), name) != loadedPluginNames.end();
}

bool LoadPlugins(const char* pluginPath, std::string parentGroup)
{
	typedef void (*FctInitPlugin)(ug::bridge::Registry*, std::string);

	bool bSuccess = true;

//	first we'll try to find all plugins in the given path
	vector<string> files;

	GetFilesInDirectory(files, pluginPath);

	bridge::Registry& reg = bridge::GetUGRegistry();

	std::string prefixStr = GetDynamicLibraryPrefix();
	std::string suffixStr = std::string(".").append(GetDynamicLibrarySuffix());

	int prefixLen = prefixStr.size();
	int suffixLen = suffixStr.size(); // includes '.'

	for(size_t i = 0; i < files.size(); ++i){

	//	extract the plugins name from the file-name
		int nameStart = prefixLen;
		int nameLength = (int)files[i].size() - suffixLen - nameStart;

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

	//	find the plugin init function
		string pluginName = files[i].substr(nameStart, nameLength);
		std::string fctName("InitUGPlugin");
		fctName.append(pluginName);

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
		loadedPlugins.push_back(loadedPluginNames);
	}

//	make sure that the registry is updated
	reg.registry_changed();

	return bSuccess;
}


bool UnloadPlugins()
{	
	// please implement
	return false;
}

}// end of namespace

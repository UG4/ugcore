/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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

bool LoadPlugins(const char* pluginPath, string parentGroup, bridge::Registry& reg, bool bPrefixGroup)
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


		DynLibHandle libHandle;
		try{
			libHandle = OpenLibrary(fullPluginName.c_str());
		}
		catch(std::string errMsg)
		{
			UG_ERR_LOG("PLUGIN-ERROR: Couldn't open plugin " << files[i] << endl);
			UG_ERR_LOG("Error Message: " << errMsg << "\n");
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

//#define DEBUG_PLUGINS
#ifdef DEBUG_PLUGINS
		UG_LOG("Loaded plugin " << pluginName << " from " << fullPluginName << "\n");
		size_t numClassesPre = reg.num_classes();
		size_t numFunctionsPre = reg.num_functions();
		UG_LOG("Call " << fctName << "... ");
#endif

		// added this for better docu generation
		// this way, we can distinguish what plugin did what classes/functions etc.
		std::string group;
		if(bPrefixGroup)
			group = std::string("(Plugin) ") + pluginName + std::string(" ") + parentGroup;
		else
			group = parentGroup;

	//	call the init method
		fctInitPlugin(&reg, group);
#ifdef DEBUG_PLUGINS
		UG_LOG("added " << reg.num_classes() - numClassesPre << " classes, " << reg.num_functions()-numFunctionsPre << " functions, "
				<< "group = " << group << "\n")
#endif

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

vector<string> GetLoadedPlugins()
{
	return loadedPluginNames;
}

}// end of namespace

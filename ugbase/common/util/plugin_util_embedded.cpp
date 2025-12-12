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
#include <vector>
#include <cstring>
#include "common/log.h"
#include "bridge/bridge.h"
#include "common/profiler/profiler.h"
#include "path_provider.h"
#include "plugin_util.h"
#include "embedded_plugins.h"	// created by cmake in the build directory

using namespace std;

namespace ug {

bool PluginLoaded(const string &name)
{
	string searchStr;
	searchStr.append(" ").append(name).append(" ");
	return (strstr(ListOfEmbeddedPlugins(), searchStr.c_str()) != nullptr);
}

bool LoadPlugins(const char*, const string &parentGroup, bridge::Registry& reg, bool bPrefixGroup)
{
	InitializeEmbeddedPlugins(&reg, parentGroup);
	return true;
}

bool UnloadPlugins()
{
	return true;
}

vector<string> GetLoadedPlugins()
{
	string plugins = string(ListOfEmbeddedPlugins());

	std::stringstream ss(plugins);
	std::istream_iterator<std::string> begin(ss);
	std::istream_iterator<std::string> end;
	std::vector<std::string> pluginNames(begin, end);
	std::copy(pluginNames.begin(), pluginNames.end(), std::ostream_iterator<std::string>(std::cout, "\n"));

	return pluginNames;
}

}// end of namespace

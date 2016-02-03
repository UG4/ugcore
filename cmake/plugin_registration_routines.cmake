# Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

################################################################################
# Register a core plugin. Registered plugins can be enabled through an option. 
# The global property ugPluginNames will contain the names of all registered plugins.
# The global property ugPluginDirs will contain the associated directories
# @param name	name of the plugin
# @param dir	directory in which the sources of the plugin lie
function(RegisterPlugin name dir)
    if (PRINT_PLUGINS_FOR_COMPLETION)
        message("-D${name}= \\")
    endif()

	# only continue if the path exists
	IF(IS_DIRECTORY ${dir})
	 	# retrieve the global property ugPluginNames and store the name in it
	    get_property(tmp GLOBAL PROPERTY ugPluginNames)
	    set_property(GLOBAL PROPERTY ugPluginNames ${tmp} "${name}")
	    
	    # same for ugPluginDirs
	    get_property(tmp GLOBAL PROPERTY ugPluginDirs)
	    set_property(GLOBAL PROPERTY ugPluginDirs ${tmp} "${dir}")

	    # add the option and set it accordingly to the defaultState parameter
	    option(${name} "Plugin: ${name}" OFF)
	    # message(STATUS "DEBUG: registered plugin ${name} at ${dir}")
	endif()
	
	mark_as_advanced(TMP_${name}_PATH)
	mark_as_advanced(TMP_${name}_PATH-NOTFOUND)	
endfunction(RegisterPlugin)
################################################################################


################################################################################
function(EnableAllPlugins)
	get_property(plugins GLOBAL PROPERTY ugPluginNames)
	foreach(plugin ${plugins})
		set(${plugin} ON CACHE BOOL "..." FORCE)
	endforeach(plugin)
endfunction(EnableAllPlugins)
################################################################################


################################################################################
function(DisableAllPlugins)
	get_property(plugins GLOBAL PROPERTY ugPluginNames)
	foreach(plugin ${plugins})
		set(${plugin} OFF CACHE BOOL "..." FORCE)
	endforeach(plugin)
endfunction(DisableAllPlugins)
################################################################################


################################################################################
function(ListPlugins)

	set(enabledPluginsStr "")
	set(disabledPluginsStr "")

	get_property(plugins GLOBAL PROPERTY ugPluginNames)
	foreach(plugin ${plugins})
		if(${${plugin}} STREQUAL "ON")
			set(enabledPluginsStr ${enabledPluginsStr} ${plugin})
		else(${${plugin}} STREQUAL "ON")
			set(disabledPluginsStr ${disabledPluginsStr} ${plugin})
		endif(${${plugin}} STREQUAL "ON")
	endforeach(plugin)
	
	message(STATUS "")
	message(STATUS "Info: Enabled plugins:")
	set(CNT 0)
	set(msg "")
	foreach(plugin ${enabledPluginsStr})
		set(msg "${msg}${plugin}, ")
		math(EXPR CNT ${CNT}+1)
		if(CNT GREATER 2)
			message(STATUS "      ${msg}")
			set(CNT 0)
			set(msg "")
		endif(CNT GREATER 2)
	endforeach(plugin ${enabledPluginsStr})
	if(CNT GREATER 0)
		message(STATUS "      ${msg}")
	endif(CNT GREATER 0)

	message(STATUS "Info: Disabled plugins:")
	set(CNT 0)
	set(msg "")
	foreach(plugin ${disabledPluginsStr})
		set(msg "${msg}${plugin}, ")
		math(EXPR CNT ${CNT}+1)
		if(CNT GREATER 2)
			message(STATUS "      ${msg}")
			set(CNT 0)
			set(msg "")
		endif(CNT GREATER 2)
	endforeach(plugin ${disabledPluginsStr})
	if(CNT GREATER 0)
		message(STATUS "      ${msg}")
	endif(CNT GREATER 0)

	message(STATUS "")
	message(STATUS "Hint: To enable/disable a plugin 'PLUGIN_NAME' use the cmake option: -DPLUGIN_NAME=ON/OFF")
	
endfunction(ListPlugins)


################################################################################
# Add subdirectories of active plugins
function(AddActivePluginSubdirectories)
	# Add subdirectories for active plugins.

	set(registeredPlugins "")
	set(registeredPluginDirs "")
	set(registeredPluginBinDirs "")
	
	get_property(plugins GLOBAL PROPERTY ugPluginNames)
	foreach(plugin ${plugins})
		set(registeredPlugins ${registeredPlugins} ${plugin})
		set(registeredPluginBinDirs ${registeredPluginBinDirs} plugins/${plugin})
	endforeach(plugin)
	
	get_property(pluginsDirs GLOBAL PROPERTY ugPluginDirs)
	foreach(pluginDir ${pluginsDirs})
		set(registeredPluginDirs ${registeredPluginDirs} ${pluginDir})
	endforeach(pluginDir)
	
	# Add directories of all enabled plugins
	if(NOT "${registeredPluginDirs}" STREQUAL "")
		list(LENGTH registeredPluginDirs numPlugins)
		math(EXPR numPlugins ${numPlugins}-1)
		
		foreach(i RANGE ${numPlugins})
			list(GET registeredPlugins ${i} pluginName)
			if(${${pluginName}} STREQUAL "ON")
				list(GET registeredPluginDirs ${i} pluginDir)
				list(GET registeredPluginBinDirs ${i} pluginBinDir)
				add_subdirectory(${pluginDir} ${pluginBinDir})
			endif(${${pluginName}} STREQUAL "ON")
		endforeach(i)
	endif()
endfunction(AddActivePluginSubdirectories)


################################################################################
function(EraseUnusedPlugins)
	get_property(plugins GLOBAL PROPERTY ugPluginNames)
	foreach(plugin ${plugins})
		if(${${plugin}} STREQUAL "OFF")
		# erase the plugin
			foreach(prefix ${CMAKE_FIND_LIBRARY_PREFIXES})
				foreach(suffix ${CMAKE_FIND_LIBRARY_SUFFIXES})
					file(REMOVE ${UG_ROOT_PATH}/bin/plugins/${prefix}${plugin}${suffix})
				endforeach(suffix) 
			endforeach(prefix) 
		endif(${${plugin}} STREQUAL "OFF")	
	endforeach(plugin)
	
endfunction(EraseUnusedPlugins)


################################################################################
function(WriteStaticPluginsHeader)

	# Collect all enabled plugins and construct initialization calls
	set(enabledPlugins "")
	set(enabledPluginsStr "")
	set(funcDeclarations "")
	set(argStr "(ug::bridge::Registry* reg, std::string parentGroup)")
	set(initCalls "")
	
	get_property(plugins GLOBAL PROPERTY ugPluginNames)
	foreach(plugin ${plugins})
		# option(${plugin} "Plugin: ${plugin}" ON)

		if(${${plugin}} STREQUAL "ON")
			set(enabledPlugins ${enabledPluginsStr} ${plugin})
			set(enabledPluginsStr "${enabledPluginsStr}${plugin} ")
			set(funcDeclarations "${funcDeclarations}extern \"C\" UG_API void InitUGPlugin_${plugin}${argStr};\n")
			set(initCalls "${initCalls}    InitUGPlugin_${plugin}(reg, parentGroup);\n")
		endif(${${plugin}} STREQUAL "ON")
	endforeach(plugin)
	
	# write a file which contains the list of active plugins
	# NOTE: the leading and closing whitespace in the string returned by InitializeEmbeddedPlugins is intentional!
	#		It allows to easily search for loaded plugins
	file(WRITE ${CMAKE_BINARY_DIR}/embedded_plugins.h
		"// THIS IS AN AUTOMATICALL GENERATED FILE. DO NOT EDIT!\n"
		"#ifndef __H__UG_STATIC_PLUGINS\n"
		"#define __H__UG_STATIC_PLUGINS\n"
		"#include <string>\n"
		"#include \"registry/registry.h\"\n\n"
		"${funcDeclarations}\n"
		"namespace ug{\n"
		"inline const char* ListOfEmbeddedPlugins()	{return \" ${enabledPluginsStr} \";}\n\n"
		"inline void InitializeEmbeddedPlugins(ug::bridge::Registry* reg, const std::string& parentGroup)\n{\n"
		"${initCalls}"
		"}\n"
		"}// end of namespace\n"
		"#endif\n")

endfunction(WriteStaticPluginsHeader)

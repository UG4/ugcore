# Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

cmake_minimum_required (VERSION 3.10...3.27.1)
message( "<ugcore>")
################################################
# ug4
################################################
project (UG4_CORE)

# set global output paths for UG4
# They are relative to the source dir for which cmake was executed.
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_SOURCE_DIR}/bin)
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_SOURCE_DIR}/bin)

include ("cmake/ug_includes.cmake")
include ("cmake/ug_cpack_config.cmake")
########################
# compile plugins
message(STATUS "${buildPlugins} buildPlugins")

if (buildPlugins)
	include ("cmake/plugin_registration_routines.cmake")
	set (pluginDir ${UG_ROOT_PATH}/plugins)
	file (GLOB children RELATIVE ${pluginDir} ${pluginDir}/*)
	foreach (child ${children})
		if (IS_DIRECTORY ${pluginDir}/${child})
			RegisterPlugin(${child} ${pluginDir}/${child})
		endif ()
	endforeach ()

	option (ENABLE_ALL_PLUGINS "If set to ON, all plugins will be enabled" OFF)
	option (DISABLE_ALL_PLUGINS "If set to ON, all plugins will be disabled" OFF)
	
	if (ENABLE_ALL_PLUGINS)
		EnableAllPlugins ()
	endif ()

	if (DISABLE_ALL_PLUGINS)
		DisableAllPlugins ()
	endif ()

	set (ENABLE_ALL_PLUGINS OFF CACHE BOOL "..." FORCE)
	set (DISABLE_ALL_PLUGINS OFF CACHE BOOL "..." FORCE)

	ListPlugins ()
	AddActivePluginSubdirectories ()

	if (buildEmbeddedPlugins)
		WriteStaticPluginsHeader ()
	else (buildEmbeddedPlugins)
		EraseUnusedPlugins ()
	endif ()
endif ()
message(STATUS "buildTools")

if (buildTools)
	set (toolDir ${UG_ROOT_PATH}/tools)
	set (FOUND_TOOLS 0)
	file (GLOB children RELATIVE ${toolDir} ${toolDir}/*)
	foreach (child ${children})
		if (IS_DIRECTORY ${toolDir}/${child})

			if (NOT ${FOUND_TOOLS})
				message(STATUS "")
				message(STATUS "Tools:")
				set (FOUND_TOOLS 1)
			endif ()

			option(${child} "  -${child}" OFF)
			
			if (${child})
				message (STATUS "  - " ${child} " (ON)" )
				add_subdirectory (${toolDir}/${child} tools/${child})
			else()
				message (STATUS "  - " ${child} " (OFF)" )
			endif ()
		endif ()
	endforeach ()
endif ()


########################
# ugbase
add_subdirectory (ugbase)

########################

if (BUILD_DOXYDOCU)
   ADD_SUBDIRECTORY (docs)
endif ()

########################
if (BUILD_ONE_LIB)
	get_property (ug4libSources GLOBAL PROPERTY ugSources)
	get_property (ug4libDependencies GLOBAL PROPERTY ugDependencies)
	get_property (ug4LinkPaths GLOBAL PROPERTY ugLinkPaths)
	get_property (ug4libIncludes GLOBAL PROPERTY ugIncludes)
	get_property (ug4libDefinitions GLOBAL PROPERTY ugDefinitions)
	get_property (ug4libLinkerFlags GLOBAL PROPERTY ugLinkerFlags)
	get_property (ug4libSingleFileCompileFlagsFiles GLOBAL PROPERTY ugSingleFileCompileFlagsFiles)
	get_property (ug4libSingleFileCompileFlagsDefs GLOBAL PROPERTY ugSingleFileCompileFlagsDefs)

	if (buildEmbeddedPlugins)
		add_definitions (${ug4libDefinitions})
		link_directories (${ug4LinkPaths})
		set (CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${ug4libLinkerFlags}")
		set (CMAKE_STATIC_LINKER_FLAGS "${CMAKE_STATIC_LINKER_FLAGS} ${ug4libLinkerFlags}")
	endif ()

	if (buildDynamicLib)
		add_definitions (-DBUILDING_DYNAMIC_LIBRARY)
		if (WIN32)
			add_definitions (-DLUA_BUILD_AS_DLL -DLUA_LIB -DLUA_CORE)
		endif ()

	   # dynamic; w/o cuda
	   add_library (${targetLibraryName} SHARED ${ug4libSources})
	   install (TARGETS ${targetLibraryName}
				LIBRARY DESTINATION lib
				COMPONENT libraries)

		#ug_add_library (${targetLibraryName} SHARED ${ug4libSources})

		if (buildCompileInfo AND NOT buildEmbeddedPlugins)
			add_dependencies (${targetLibraryName} CompileInfo)
			target_link_libraries (${targetLibraryName} CompileInfo)
		endif ()

	else()

		add_library (${targetLibraryName} ${ug4libSources})

		if (buildCompileInfo AND NOT buildEmbeddedPlugins)
			add_dependencies(${targetLibraryName} CompileInfo_s)
			target_link_libraries(${targetLibraryName} CompileInfo_s)

		endif ()
	endif ()

  target_include_directories (${targetLibraryName} PUBLIC ${ug4libIncludes})
  target_link_libraries (${targetLibraryName} ${ug4libDependencies} ${linkLibraries})

  target_enable_sanitizers (${targetLibraryName} ${USE_SANITIZER})

	if (buildEmbeddedPlugins)
		# if we build with embedded plugins, the compile-info is also embedded
		# into the single big target lib. We thus have to add a dependency to
		# the custom updateCompileInfo target declared in compile_info/CMakeLists.txt.
    if (UNIX AND buildCompileInfo)
		  add_dependencies(${targetLibraryName} updateCompileInfo)
    endif ()

		# apply single file definitions
		list (LENGTH ug4libSingleFileCompileFlagsFiles len)
		if (${len} GREATER 0)
			math (EXPR len "${len} - 1")
			foreach (idx RANGE ${len})
				list (GET ug4libSingleFileCompileFlagsFiles ${idx} file)
				list (GET ug4libSingleFileCompileFlagsDefs ${idx} flags)
				set_property (SOURCE ${file} APPEND PROPERTY COMPILE_FLAGS "${flags}")
				#message ("${file}: ${flags}")
			endforeach ()
		endif ()
	endif ()

endif ()
message( "</ugcore>")

get_property(_t DIRECTORY PROPERTY BUILDSYSTEM_TARGETS)
message(STATUS "Directory targets: ${_t}")
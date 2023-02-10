# Copyright (c) 2012-2013:  G-CSC, Goethe University Frankfurt
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

cmake_minimum_required(VERSION 2.8.12)

################################################
# ug4 - compile-info
################################################

project(P_UG4_COMPILE_INFO)

# set global output paths for UG4
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("${UG_ROOT_CMAKE_PATH}/ug_includes.cmake")

#	Name of the file that will contain the compile-info variables
set(infoFileName ${CMAKE_BINARY_DIR}/compile_info_vars.cpp)

#	we write this file to make sure that it exists...
#	Note that on UNIX systems the file is updated with each build through
#	the add_custom_build_target and add_custom_command calls below. 
#	on windows we currently don't support updated compile-infos.
#	check out, whether this can be improved using special windows-scripting...

#set(infoFileFound "infoFileFound-NOTFOUND")
#find_file(infoFileFound compile_info_vars.cpp PATHS ${CMAKE_BINARY_DIR})

#if(${infoFileFound} STREQUAL "infoFileFound-NOTFOUND")
	file(WRITE ${infoFileName} 
		"namespace ug\n"
		"{\n"
		"const char *UG_SVN_REVISION=\"UNKNOWN\";\n"
		"const char *UG_GIT_REVISION=\"UNKNOWN\";\n"
		"const char *UG_BUILD_HOST=\"UNKNOWN\";\n"
		"const char *UG_COMPILE_DATE=\"UNKNOWN\";\n"
		"}\n")
#endif(${infoFileFound} STREQUAL "infoFileFound-NOTFOUND")	

if(UNIX AND buildCompileInfo)
#	The custom target is always considered out of date and will be rebuilt
#	every time make is executed.
	add_custom_target(updateCompileInfo)
#	The build of the custom target triggers this custom command, which updates
#	the compile_info.cpp file.
	add_custom_command(TARGET updateCompileInfo PRE_BUILD COMMAND bash ${UG_ROOT_PATH}/ugcore/scripts/shell/generate_compile_info ${UG_ROOT_PATH} ${infoFileName})
endif(UNIX AND buildCompileInfo)

#	Export sources to the file which contains the info-access-methods.
#	Note that this is not the auto-generated file which contains
#	the actual variable definitions.
EXPORTSOURCES(ugbase/compile_info compile_info.cpp)

if(buildEmbeddedPlugins OR NOT buildCompileInfo)
	EXPORTSOURCES("" ${infoFileName})
#	Note that a dependency to updateCompileInfo is created in the main CMakeLists.txt
#	file for this case.

else(buildEmbeddedPlugins OR NOT buildCompileInfo)
	if(buildDynamicLib)
		add_library(CompileInfo SHARED ${infoFileName})
		install(TARGETS CompileInfo
			ARCHIVE
  			DESTINATION lib)
		set(compileInfoLibName CompileInfo)
	else(buildDynamicLib)
	   # no dynamic lib
		add_library(CompileInfo_s STATIC ${infoFileName})
		set(compileInfoLibName CompileInfo_s)
	endif(buildDynamicLib)
	
	if(UNIX AND buildCompileInfo)
	#	The dependency makes sure, that the compile_info file is written before
	# 	the CompileInfo lib is being rebuilt.
		add_dependencies(${compileInfoLibName} updateCompileInfo)
	endif(UNIX AND buildCompileInfo)
endif(buildEmbeddedPlugins OR NOT buildCompileInfo)

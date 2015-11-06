# Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

# This file is included e.g. by ug4/trunk/CMakeLists.txt and defines the
# content of the cmake-root-file.

cmake_minimum_required(VERSION 2.6)

################################################
# ug4
################################################
project(P_UG4)

# set global output paths for UG4
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("cmake/ug_includes.cmake")

########################
# compile plugins
if(buildPlugins)
   add_subdirectory(plugins)
endif(buildPlugins)

########################
# ugbase
add_subdirectory(ugbase)

########################
if(BUILD_UGDOCU)
    add_subdirectory(apps/ugdocu)
endif(BUILD_UGDOCU)

IF(BUILD_DOXYDOCU)
   ADD_SUBDIRECTORY(docs)
ENDIF(BUILD_DOXYDOCU)

########################
if(BUILD_ONE_LIB)
	get_property(ug4libSources GLOBAL PROPERTY ugSources)
	get_property(ug4libDependencies GLOBAL PROPERTY ugDependencies)
	get_property(ug4libIncludes GLOBAL PROPERTY ugIncludes)
	get_property(ug4libDefinitions GLOBAL PROPERTY ugDefinitions)
	
	if(buildEmbeddedPlugins)
   add_definitions(${ug4libDefinitions})
   include_directories(${ug4libIncludes})
		link_libraries(${ug4libDependencies})
	endif(buildEmbeddedPlugins)
	
	if(buildDynamicLib)
		add_definitions(-DBUILDING_DYNAMIC_LIBRARY)
		if(WIN32)
			add_definitions(-DLUA_BUILD_AS_DLL -DLUA_LIB -DLUA_CORE)
		endif(WIN32)
		if(CUDA_FOUND)
		   CUDA_add_library(${targetLibraryName} SHARED ${ug4libSources})
		else(CUDA_FOUND)
		   add_library(${targetLibraryName} SHARED ${ug4libSources})
		endif(CUDA_FOUND)
		#ug_add_library(${targetLibraryName} SHARED ${ug4libSources})
		
		if(buildCompileInfo AND NOT buildEmbeddedPlugins)
			add_dependencies(${targetLibraryName} CompileInfo)
			target_link_libraries(${targetLibraryName} CompileInfo)
		endif(buildCompileInfo AND NOT buildEmbeddedPlugins)

   else(buildDynamicLib)
       if(CUDA_FOUND)
           CUDA_add_library(${targetLibraryName} ${ug4libSources})
       else(CUDA_FOUND)
           add_library(${targetLibraryName} ${ug4libSources})
        endif(CUDA_FOUND)
       
      
      if(buildCompileInfo AND NOT buildEmbeddedPlugins)
         add_dependencies(${targetLibraryName} CompileInfo_s)
         target_link_libraries(${targetLibraryName} CompileInfo_s)
      endif(buildCompileInfo AND NOT buildEmbeddedPlugins)
      
   endif(buildDynamicLib)
   
   if(buildEmbeddedPlugins)
   #  if we build with embedded plugins, the compile-info is also embedded
   #  into the single big target lib. We thus have to add a dependency to
   #  the custom updateCompileInfo target declared in compile_info/CMakeLists.txt.
      add_dependencies(${targetLibraryName} updateCompileInfo)
   endif(buildEmbeddedPlugins)
   
endif(BUILD_ONE_LIB)

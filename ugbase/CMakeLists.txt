################################################
# created by Sebastian Reiter, Andreas Vogel
# s.b.reiter@googlemail.com
#
# ugbase
################################################

cmake_minimum_required(VERSION 2.6)

project(P_UGBASE)

# set global output paths for UGBASE
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(Sources	ug.cpp)

if(BUILD_ONE_LIB)
	EXPORTSOURCES(. ${Sources})
else(BUILD_ONE_LIB)
	add_library(ugbase ${Sources})
endif(BUILD_ONE_LIB)


########################
# common
add_subdirectory(common)

########################
# node_tree
add_subdirectory(node_tree)

########################
# lib_grid
add_subdirectory(lib_grid)

########################
# lib_registry
add_subdirectory(registry)

########################
# lib_ugbridge
add_subdirectory(ug_bridge)

########################
# lib_script
add_subdirectory(ug_script)

########################
# lib_algebra and lib_discretization
if(buildAlgebra)
	add_subdirectory(lib_algebra)	
	add_subdirectory(lib_discretization)
endif(buildAlgebra)


########################
# pcl
if(PARALLEL)
	add_subdirectory(pcl)
endif(PARALLEL)


########################
# vrl_bindings
if(buildForVRL)
    add_subdirectory(bindings_vrl)
endif(buildForVRL)

########################
# ug4 library
if(BUILD_ONE_LIB)
	get_property(ug4libSources GLOBAL PROPERTY ugSources)
	if(buildDynamicLib)
		add_definitions(-DBUILDING_DYNAMIC_LIBRARY)
		if(WIN32)
			add_definitions(-DLUA_BUILD_AS_DLL -DLUA_LIB -DLUA_CORE)
		endif(WIN32)
		add_library(${targetLibraryName} SHARED ${ug4libSources})
	else(buildDynamicLib)
		add_library(${targetLibraryName} ${ug4libSources})
	endif(buildDynamicLib)
endif(BUILD_ONE_LIB)

########################
# ug4 shell
if(buildUGShell)
    add_subdirectory(ug_shell)
endif(buildUGShell)


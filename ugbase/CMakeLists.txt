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
# lib_grid
if(buildGrid)
	add_subdirectory(lib_grid)
endif(buildGrid)

########################
# lib_registry
if(buildRegistry)
	add_subdirectory(registry)
endif(buildRegistry)

########################
# bridge
if(buildBridge)
	add_subdirectory(bridge)
endif(buildBridge)

########################
# bindings
if(buildBindings)
	add_subdirectory(bindings/lua)
	if(buildForVRL)
    	add_subdirectory(bindings/vrl)
	endif(buildForVRL)
endif(buildBindings)

########################
# lib_algebra
if(buildAlgebra)
	add_subdirectory(lib_algebra)	
endif(buildAlgebra)

########################
# lib_disc
if(buildDisc)
	add_subdirectory(lib_disc)
endif(buildDisc)


########################
# pcl
if(PARALLEL)
	add_subdirectory(pcl)
endif(PARALLEL)

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
#if(buildUGShell)
#    add_subdirectory(ug_shell)
#endif(buildUGShell)


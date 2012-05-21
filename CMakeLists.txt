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

include("ug_cmake_includes.txt")



########################
# ugbase

add_subdirectory(ugbase)

# compile core plugins
if(buildCorePlugins)
	add_subdirectory(plugins)
endif(buildCorePlugins)

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

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
# compile plugins
if(buildPluginSystem)
	add_subdirectory(plugins)
endif(buildPluginSystem)


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
	
	if(buildDynamicLib)
		add_definitions(-DBUILDING_DYNAMIC_LIBRARY)
		if(WIN32)
			add_definitions(-DLUA_BUILD_AS_DLL -DLUA_LIB -DLUA_CORE)
		endif(WIN32)
		add_library(${targetLibraryName} SHARED ${ug4libSources})

		if(buildCompileInfo AND NOT buildEmbeddedPlugins)
			add_dependencies(${targetLibraryName} CompileInfo)
			target_link_libraries(${targetLibraryName} CompileInfo)
		endif(buildCompileInfo AND NOT buildEmbeddedPlugins)

	else(buildDynamicLib)
		add_library(${targetLibraryName} ${ug4libSources})
		
		if(buildCompileInfo AND NOT buildEmbeddedPlugins)
			add_dependencies(${targetLibraryName} CompileInfo_s)
			target_link_libraries(${targetLibraryName} CompileInfo_s)
		endif(buildCompileInfo AND NOT buildEmbeddedPlugins)
		
	endif(buildDynamicLib)
	
	if(buildEmbeddedPlugins)
	#	if we build with embedded plugins, the compile-info is also embedded
	#	into the single big target lib. We thus have to add a dependency to
	#	the custom updateCompileInfo target declared in compile_info/CMakeLists.txt.
		add_dependencies(${targetLibraryName} updateCompileInfo)
	endif(buildEmbeddedPlugins)
	
endif(BUILD_ONE_LIB)




cmake_minimum_required(VERSION 2.6)

################################################
# ug4 - compile-info
################################################

project(P_UG4_COMPILE_INFO)

# set global output paths for UG4
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("../../cmake/ug_includes.cmake")

#	Name of the file that will contain the compile-info variables
set(infoFileName ${CMAKE_BINARY_DIR}/compile_info_vars.cpp)

#	we write this file to make sure that it exists...
#	Note that on UNIX systems the file is updated with each build through
#	the add_custom_build_target and add_custom_command calls below. 
#	on windows we currently don't support updated compile-infos.
#	check out, whether this can be improved using special windows-scripting...

unset(infoFileFound CACHE)
find_file(infoFileFound compile_info_vars.cpp PATHS ${CMAKE_BINARY_DIR})

if(${infoFileFound} STREQUAL "infoFileFound-NOTFOUND")
	file(WRITE ${infoFileName} 
		"namespace ug\n"
		"{\n"
		"const char *UG_SVN_REVISION=\"UNKNOWN\";\n"
		"const char *UG_BUILD_HOST=\"UNKNOWN\";\n"
		"const char *UG_COMPILE_DATE=\"UNKNOWN\";\n"
		"}\n")
endif(${infoFileFound} STREQUAL "infoFileFound-NOTFOUND")	

if(UNIX AND buildCompileInfo)
#	The custom target is always considered out of date and will be rebuilt
#	every time make is executed.
	add_custom_target(updateCompileInfo)
#	The build of the custom target triggers this custom command, which updates
#	the compile_info.cpp file.
	add_custom_command(TARGET updateCompileInfo PRE_BUILD COMMAND ${UG_ROOT_PATH}/scripts/shell/generate_compile_info ${UG_ROOT_PATH} ${infoFileName})
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
		set(compileInfoLibName CompileInfo)
	else(buildDynamicLib)
		add_library(CompileInfo_s STATIC ${infoFileName})
		set(compileInfoLibName CompileInfo_s)
	endif(buildDynamicLib)
	
	if(UNIX AND buildCompileInfo)
	#	The dependency makes sure, that the compile_info file is written before
	# 	the CompileInfo lib is being rebuilt.
		add_dependencies(${compileInfoLibName} updateCompileInfo)
	endif(UNIX AND buildCompileInfo)
endif(buildEmbeddedPlugins OR NOT buildCompileInfo)

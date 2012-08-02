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

include("../../ug_cmake_includes.txt")

#	we write this file to make sure that it exists...
#	Note that on UNIX systems the file is updated with each build through
#	the add_custom_build_target and add_custom_command calls below. 
#	on windows we currently don't support updated compile-infos.
#	check out, whether this can be improved using special windows-scripting...
file(WRITE ${CMAKE_BINARY_DIR}/compile_info.cpp 
	"namespace ug\n"
	"{\n"
	"const char *SVN_REVISION=\"UNKNOWN\";\n"
	"const char *BUILD_HOST="\"${CMAKE_SYSTEM_NAME}\";\n"
	"const char *COMPILE_DATE="\"UNKNOWN\";\n"
	"}\n")

IF(UNIX)
#	The custom target is always considered out of date and will be rebuilt
#	every time make is executed.
	add_custom_target(updateCompileInfo)
#	The build of the custom target triggers this custom command, which updates
#	the compile_info.cpp file.
	add_custom_command(TARGET updateCompileInfo PRE_BUILD COMMAND ${UG_ROOT_PATH}/scripts/shell/generate_compile_info ${UG_ROOT_PATH} ${CMAKE_BINARY_DIR}/compile_info.cpp)
ENDIF(UNIX)

if(buildEmbeddedPlugins)
	EXPORTSOURCES(/ ${CMAKE_BINARY_DIR}/compile_info.cpp)
#	Note that a dependency to updateCompileInfo is created in the main CMakeLists.txt
#	file for this case.

else(buildEmbeddedPlugins)
#	this clearly is a dirty hack...
#	we have to prevent that the library tries to link against itself.
#	The benefit of setting the linker flags in the ug_cmake_includes file however
#	currently outweights this ugly stuff...
#	Note that this has no influence on other libraries...
	set(CMAKE_SHARED_LINKER_FLAGS)
	set(CMAKE_EXE_LINKER_FLAGS)
	set(CMAKE_MODULE_LINKER_FLAGS)
		
	if(buildDynamicLib)
		add_library(CompileInfo SHARED ${CMAKE_BINARY_DIR}/compile_info.cpp)
		set(compileInfoLibName CompileInfo)
	else(buildDynamicLib)
		add_library(CompileInfo_s STATIC ${CMAKE_BINARY_DIR}/compile_info.cpp)
		set(compileInfoLibName CompileInfo_s)
	endif(buildDynamicLib)
	
	IF(UNIX)
	#	The dependency makes sure, that the compile_info file is written before
	# 	the CompileInfo lib is being rebuilt.
		add_dependencies(${compileInfoLibName} updateCompileInfo)
	ENDIF(UNIX)
endif(buildEmbeddedPlugins)

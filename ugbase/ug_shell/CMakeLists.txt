cmake_minimum_required(VERSION 2.6)

####
# ugscript Library
####

project(P_UGSCRIPT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("../../ug_cmake_includes.txt")

set(shelltypeOptions "linenoise, readline, simple")
set(shelltypeDefault "linenoise")

# In order to avoid build errors, we'll perform some more checks, to determine
# the best matching shelltype for the given platform.
#todo: This can be further elaborated...
if(NOT UNIX)
	# set shelltypeDefault to "simple". linenoise most likely wouldn't work, too.
	set(shelltypeDefault "simple")
endif(NOT UNIX)


# this option is a pseudo-option, since it uses a string instead of (ON / OFF)
option(SHELLTYPE "Set the type of the shell. Valid options are: "${targetOptions})

if(NOT SHELLTYPE)
	set(SHELLTYPE ${shelltypeDefault})
endif(NOT SHELLTYPE)

# output the current shell-type
message(STATUS "")
message(STATUS "Info: SHELLTYPE: "${SHELLTYPE}" (options are: "${shelltypeOptions}")")


# set the main sources
set(srcUGShell	    ugshell_main.cpp
                    completion.cpp)
                    #signal_handler.cpp)


if("${SHELLTYPE}" STREQUAL "linenoise")
	add_definitions(-DUG_USE_LINENOISE)
	set(srcUGShell	${srcUGShell} externals/linenoise/linenoise.cpp)
	
elseif("${SHELLTYPE}" STREQUAL "readline")
	find_library(libReadline readline)
	if(libReadline)
		add_definitions(-DUG_USE_READLINE)
		message(STATUS "Info: Readline is GPL licensed. Do not redistribute the resulting application!")
	else(libReadline)
		message(FATAL_ERROR "Readline could not be found. Aborting.")
	endif(libReadline)

elseif("${SHELLTYPE}" STREQUAL "simple")
	# Nothing to do here...
	
else("${SHELLTYPE}" STREQUAL "linenoise")
	message(FATAL_ERROR "Invalid SHELLTYPE specified. Aborting.")
endif("${SHELLTYPE}" STREQUAL "linenoise")


message(STATUS)

remove_definitions(-DBUILDING_DYNAMIC_LIBRARY)
if(buildDynamicLibrary)
	add_definitions(-DIMPORT_DYNAMIC_LIBRARY)
endif(buildDynamicLibrary)

remove_definitions(-DLUA_BUILD_AS_DLL)

add_executable(${targetExecutableName} ${srcUGShell})

if("${SHELLTYPE}" STREQUAL "readline")
	target_link_libraries(${targetExecutableName} readline)
endif("${SHELLTYPE}" STREQUAL "readline")

if(STATIC)
# if we're compiling with the new cmake_includes, then we have to link
# against the ug4_s static library
	target_link_libraries(${targetExecutableName} ${targetLibraryName})
else(STATIC)
# link the dynamic library
	target_link_libraries (${targetExecutableName} ${targetLibraryName})
endif(STATIC)

if(BUILD_UGDOCU)
	add_custom_target(updateUGDocu ALL )
	add_custom_command(TARGET updateUGDocu
						 PRE_BUILD
						 COMMAND "${UG_ROOT_PATH}/bin/ugdocu" -silent
						 WORKING_DIRECTORY ${UG_ROOT_PATH}/bin)
	add_dependencies(updateUGDocu ${TARGET} ugdocu)
endif(BUILD_UGDOCU)
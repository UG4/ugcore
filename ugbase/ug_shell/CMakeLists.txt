cmake_minimum_required(VERSION 2.6)

####
# ugscript Library
####

project(P_UGSCRIPT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

# include("../../ug_cmake_includes.txt")

# print the different options.
message(STATUS)
message(STATUS "Please note that several options exist:")
message(STATUS "  -DNO_LINENOISE=OFF: BSD-licensed linenoise is used (this is the default).")
message(STATUS "                      Adds support for a history, and path-completion.")
message(STATUS "                      Additionally type and member completion is featured.")
message(STATUS "  -DUSE_READLINE=ON: GPL-licensed readline is used.")
message(STATUS "                     Adds support for a history, and path-completion.")
message(STATUS "                     Note that you may only use this option for private applications.")
message(STATUS "                     You may not redistribute applications that use this option.")
message(STATUS "")
message(STATUS "Note that those libraries do not work for all platforms")
message(STATUS)

set(srcUGShell	    src/main.cpp
                    src/completion.cpp)

# find the readline libraries
if(USE_READLINE)
	message(STATUS "using readline. DON'T REDISTRIBUTE THE RESULTING APPLICATION!")
	find_library(libReadline readline)
	if(libReadline)
		add_definitions(-DUG_USE_READLINE)
	endif(libReadline)
else(USE_READLINE)
	if(NOT NO_LINENOISE)
		message(STATUS "INFO: Using linenoise. DON'T FORGET TO MENTION THIS IN THE DOCUMENTATION.")
		add_definitions(-DUG_USE_LINENOISE)
		set(srcUGShell	${srcUGShell} src/externals/linenoise/linenoise.cpp)
	endif(NOT NO_LINENOISE)
endif(USE_READLINE)
message(STATUS)

remove_definitions(-DBUILDING_DYNAMIC_LIBRARY)
remove_definitions(-DLUA_BUILD_AS_DLL)
add_executable(${targetExecutableName} ${srcUGShell})

if(libReadline)
	target_link_libraries(${targetExecutableName} readline)
endif(libReadline)

if(STATIC)
# if we're compiling with the new cmake_includes, then we have to link
# against the ug4_s static library
	if(USE_NEW_CMAKE_INCLUDES)
		target_link_libraries(${targetExecutableName} ${targetLibraryName})
	endif(USE_NEW_CMAKE_INCLUDES)
else(STATIC)
# link the dynamic library
	target_link_libraries (${targetExecutableName} ${targetLibraryName})
endif(STATIC)

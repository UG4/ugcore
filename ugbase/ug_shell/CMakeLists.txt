# Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

cmake_minimum_required(VERSION 2.6)

####
# ugscript Library
####

project(P_UGSCRIPT)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include("../../cmake/ug_includes.cmake")

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
option(SHELLTYPE "Set the type of the shell. Valid options are: ${targetOptions}")

if(NOT SHELLTYPE)
	set(SHELLTYPE ${shelltypeDefault})
endif(NOT SHELLTYPE)

# output the current shell-type
message(STATUS "")
message(STATUS "Info: SHELLTYPE: ${SHELLTYPE} (options are: ${shelltypeOptions})")


# set the main sources
set(srcUGShell	    ugshell_main.cpp
                    shell.cpp)
                    #algebra_test.cpp)
                    #signal_handler.cpp)

if(CLOCK_FIX)
    set(srcUGShell	    ${srcUGShell} clock_fix.cpp)
    message(STATUS "Info: USING THE CLOCK_FIX WORKAROUND (fixes undefined reference to clock_gettime and clock_getres, see ugbase/ug_shell/clock_fix.cpp).")
endif(CLOCK_FIX)

if("${SHELLTYPE}" STREQUAL "linenoise")
	add_definitions(-DUG_USE_LINENOISE)
	if(POSIX)
		set(srcUGShell	    ${srcUGShell} completion.cpp)
	endif(POSIX)
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

if(STATIC_BUILD)
    # some compiler introduce -Wl,-Bdynamic and -Wl,-Bstatic in order to 
    # distinguish linkage type for libraries. In cmake, dynamic is default.
    # Setting the options below to ON forces cmake to start and end with the 
    # static linker flag when specifying library locations. This is needed, 
    # e.g. when a compiler wrapper (like scalasca) adds its own libraries at the
    # end and the libraries are available  as shared and (!) static on the 
    # system. In that case, without the oprtions below, the dynamic lib would 
    # be selected. Setting options to ON lets cmake choose the static lib.
    set_target_properties(${targetExecutableName} PROPERTIES LINK_SEARCH_START_STATIC ON)
    set_target_properties(${targetExecutableName} PROPERTIES LINK_SEARCH_END_STATIC ON)
    
    # we have to link against the ug4_s static library
	target_link_libraries(${targetExecutableName} ${targetLibraryName})
else(STATIC_BUILD)
    # link the dynamic library
	target_link_libraries (${targetExecutableName} ${targetLibraryName})
endif(STATIC_BUILD)

include(../../cmake/build_ugdocu.cmake)

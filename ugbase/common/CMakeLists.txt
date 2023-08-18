# Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

cmake_minimum_required(VERSION 2.8.12...3.27.1)

################################################
# common
################################################

project(P_COMMON)

include("../../cmake/ug_includes.cmake")

set(srcCommon	log.cpp
				debug_id.cpp
				assert.cpp
				error.cpp
				serialization.cpp
				progress.cpp
				allocators/small_object_allocator.cpp
				util/base64_file_writer.cpp
				util/binary_buffer.cpp
				util/binary_stream.cpp
				util/demangle.cpp
				util/crc32.cpp
        		util/file_util.cpp
        		util/loader/loader_util.cpp
				util/loader/loader_obj.cpp
				util/message_hub.cpp
				util/ostream_buffer_splitter.cpp
				util/parameter_parsing.cpp
        		util/string_util.cpp
				util/variant.cpp
				util/histogramm.cpp
				util/number_util.cpp
				math/math_vector_matrix/math_matrix.cpp
				math/math_vector_matrix/math_vector.cpp
				math/misc/tri_box.cpp
				math/misc/tritri.cpp
				math/misc/lineintersect_utils.cpp
				math/misc/eigenvalues.cpp
				math/misc/math_util.cpp
				math/misc/orthopoly.cpp)
				
if(PROFILE_MEMORY)
    message(STATUS "Info: Using Memory Profiler (disable with -DPROFILE_MEMORY=OFF).")
    set(srcCommon ${srcCommon} profiler/memtracker.cpp)
else(PROFILE_MEMORY)
    set(srcCommon ${srcCommon} profiler/memtracker_disabled.cpp)
endif(PROFILE_MEMORY)

set(sources ${srcCommon})

if(buildPlugins)
	if(buildEmbeddedPlugins)
		set(sources	${sources} util/plugin_util_embedded.cpp)
	else(buildEmbeddedPlugins)
	    set(sources	${sources} util/plugin_util_dynamic.cpp)
	endif(buildEmbeddedPlugins)
endif(buildPlugins)

if(UG_PROFILER)
	set(srcProfiler	profiler/profilenode_management.cpp
	                profiler/runtime_profile_info.cpp)	
	set(sources ${sources} ${srcProfiler})
endif(UG_PROFILER)

if(UG_PROFILER_SHINY)
	set(srcShiny	profiler/src/ShinyManager.cpp
					profiler/src/ShinyNode.cpp
					profiler/src/ShinyNodePool.cpp
					profiler/src/ShinyOutput.cpp
					profiler/src/ShinyTools.cpp)	
	set(sources ${sources} ${srcShiny})
endif(UG_PROFILER_SHINY)

if(UG_CPU_FREQ)
	set(freqShiny	profiler/freq_adapt.cpp)
	set(sources ${sources} ${freqShiny})    
endif(UG_CPU_FREQ)

if(SHINY_CALL_LOGGING)
    set(sources ${sources} profiler/shiny_call_logging.cpp)
endif(SHINY_CALL_LOGGING)

# add support for UGProfileNode any case
set(sources ${sources} profiler/profile_node.cpp)

################################################################################
# Platform dependend code
################################################################################
if(UNIX)
	set(sources ${sources}	util/os_dependent_impl/file_util_posix.cpp)
	
	#if(NOT STATIC)
		set(sources ${sources}	util/os_dependent_impl/dynamic_library_util_unix.cpp)
	#endif(NOT STATIC)
	
	if(APPLE)
		set(sources ${sources} util/os_dependent_impl/os_info_osx.cpp)
	else(APPLE)
		set(sources ${sources} util/os_dependent_impl/os_info_linux.cpp)
	endif(APPLE)

elseif(WIN32)
	set(sources ${sources}	util/os_dependent_impl/os_info_win.cpp)
	
	#if(NOT STATIC)
		set(sources ${sources}	util/os_dependent_impl/dynamic_library_util_win.cpp)
	#endif(NOT STATIC)
	
	#if(MINGW)
	#	set(sources ${sources} util/os_dependent_impl/file_util_posix.cpp)
	#else(MINGW)
		set(sources ${sources} util/os_dependent_impl/file_util_win.cpp)
	#endif(MINGW)
	
else(UNIX)
	message(STATUS "YOUR OS MAY NOT BE FULLY SUPPORTED (NOT UNIX???). File functions may be not working.")
	set(sources ${sources}	util/os_dependent_impl/file_util_posix.cpp)
	set(sources ${sources} util/os_dependent_impl/os_info_linux.cpp)
endif(UNIX)


if (APPLE)
	set(sources ${sources} util/mem_info.cpp util/mem_info_apple.cpp)
endif (APPLE)

if (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")
	set(sources ${sources} util/mem_info.cpp util/mem_info_linux.cpp)
endif (${CMAKE_SYSTEM_NAME} STREQUAL "Linux")



################################################################################
# node_tree
################################################################################
add_subdirectory(node_tree)


################################################################################
#	add sources to library
################################################################################
if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/common ${sources})
else(BUILD_ONE_LIB)
	add_library(common ${sources})
endif(BUILD_ONE_LIB)

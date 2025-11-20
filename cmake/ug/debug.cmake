# Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
# Authors: Sebastian Reiter, Martin Rupp
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

# included from ug_includes.cmake
if (DEBUG)
	add_definitions (-DUG_DEBUG)
	# if no build type set, add default debug flags
	if (NOT CMAKE_BUILD_TYPE OR "${CMAKE_BUILD_TYPE}" STREQUAL "None" OR
		"${CMAKE_BUILD_TYPE}" STREQUAL "")
		# use cmake standard flags for debug builds (-g)
		add_cxx_flag (${CMAKE_CXX_FLAGS_DEBUG})
		# if user specified additional debug flags add them
		if (DEBUG_FORMAT)
			add_cxx_flag (${DEBUG_FORMAT})
			message (STATUS "Info: Debug Information is ${DEBUG_FORMAT}") 
		endif (DEBUG_FORMAT)
	elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Release")
		message (WARNING "CMAKE_BUILD_TYPE type set to Release, but DEBUG build wanted.")
	elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "MinSizeRel")
		message (WARNING "CMAKE_BUILD_TYPE type set to MinSizeRel, but DEBUG build wanted.")
	endif ()
	
	# This code would enable strict bounds checking for STL objects like in vector::operator [] .
	# however, GLIBCXX_DEBUG and strstream don't work together on mac (bug: http://bit.ly/cH78bC). 
	# when this bug is fixed, one could set those flags (or similar) depending on the compiler.
	
	# I disabled the following definitions, since they cause several additional problems.
	# First they also cause problems on MinGW builds.
	# Secondly they have to be defined also for all executables, which link against
	# ug. Since they are well hidden, this may cause quite some headache.
	# One should think about a flag DEBUG_STL or something like this, which explicitly
	# enables those flags. (sreiter)
	#if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" AND NOT APPLE)
	#	add_definitions (-D_GLIBCXX_DEBUG=1 -D_GLIBCXX_DEBUG_PEDANTIC=1)
	#endif ()
	
	# Since they are still useful I added a manual flag to enable them.
	# Note that this only works if DEBUG=ON.
	if (DEBUG_STL)
		add_definitions (-D_GLIBCXX_DEBUG=1 -D_GLIBCXX_DEBUG_PEDANTIC=1)
		message (STATUS "Info: Debugging STL")
	endif ()

	if (PARALLEL)
           add_definitions (-DLG_DISTRIBUTION_DEBUG)
	endif ()	
	
else ()
	# add release definitions
	add_definitions (-DBOOST_UBLAS_NDEBUG)
	if (NOT CMAKE_BUILD_TYPE)
		# if DEBUG=OFF also add standard cmake release flags (-O3 -DNDEBUG)
		add_cxx_flag (${CMAKE_CXX_FLAGS_RELEASE})
	elseif ("${CMAKE_BUILD_TYPE}" STREQUAL "Debug")
		message (WARNING "CMAKE_BUILD_TYPE type is set to Debug, but DEBUG=OFF. Leads to strange cflags!")
	endif ()

	# compiler specific release flags
	if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
		add_cxx_flag ("-funroll-loops -ftree-vectorize")
	elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Cray")
		add_cxx_flag ("-hipa5 -hunroll2")
	endif ()

endif ()

if (DEBUG_LOGS)
	add_definitions (-DUG_ENABLE_DEBUG_LOGS)
endif ()

if (UG_INTERSECTION_DEBUG_PATH)
	add_definitions (-DUG_INTERSECTION_DEBUG_PATH=\"${UG_INTERSECTION_DEBUG_PATH}\")
endif ()

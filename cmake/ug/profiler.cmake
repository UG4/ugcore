# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Martin Rupp
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
########################################
# PROFILER

if( NOT "${PROFILER}" STREQUAL "Shiny" AND SHINY_CALL_LOGGING)
    message(FATAL_ERROR " Shiny Call Logging activated but not Shiny. Use cmake -DPROFILER=Shiny ..")
endif( NOT "${PROFILER}" STREQUAL "Shiny" AND SHINY_CALL_LOGGING)


if(NOT "${PROFILER}" STREQUAL "None")
    if("${PROFILER}" STREQUAL "Shiny")
    	add_definitions(-DUG_PROFILER_SHINY)    
     	set(UG_PROFILER_SHINY ON)               # add Cmake variable
     	   
     	# PROFILE_BRIDGE
        if(SHINY_CALL_LOGGING)
        	add_definitions(-DSHINY_CALL_LOGGING)
        	message(" -- Info: Shiny Call Logging activated.")
        endif(SHINY_CALL_LOGGING)
             	
        
    # Scalasca
    elseif("${PROFILER}" STREQUAL "Scalasca")
        find_package(Scalasca)
        if(SCALASCA_FOUND)
            message("-- Info: Scalasca: using scalasca command: ${SCALASCA_COMMAND}")
            message("-- Info: Scalasca: using include dir: ${SCALASCA_INCLUDE_DIR}")
            message("-- Info: Scalasca: using cxx flags: ${SCALASCA_USER_CFLAGS}")
            message("-- Info: Scalasca: check that compiler prefix is set!")
            message("--       Use: CC=\"scalasca -instrument -comp=none -user mpicc\" CXX=\"scalasca -instrument -comp=none -user mpicxx\" cmake ...")
            message("--       If not used: remove build completely and rerun cmake with prefix.")
        else(SCALASCA_FOUND)
        	message(FATAL_ERROR "PROFILER: ${PROFILER}: Cannot find required "
        	        "binary scalasca/kconfig. Make sure that PATH contains "
        	        "scalasca and kconfig executable.")
        endif(SCALASCA_FOUND)

        # we assume the user to build using:
        # CXX="scalasca -instrument -comp=none -user c++"
        # It would be desirable to set the compiler here, but somehow I did not
        # find a solution to that issue. Maybe we could include this file before
        # the project() command to solve the issue, but I don't know what effects
        # that may have.
        # I also tried:

        # a) write compiler wrapper - cannot be set afterwards
#       	file(WRITE ${CMAKE_BINARY_DIR}/scalasca_mpicxx
#    		"#! /bin/sh"
#    		"# THIS IS AN AUTOMATICALL GENERATED FILE. DO NOT EDIT!\n"
#    		"${SCALASCA_COMMAND} -instrument -comp=none -user ${CMAKE_CXX_COMPILER}")
#        set(CMAKE_CXX_COMPILER "${CMAKE_BINARY_DIR}/scalasca_mpicxx")

         # b) resetting compiler: does not work, only allowed before PROJECT() command
#        if(NOT ("${CMAKE_CXX_COMPILER}" STREQUAL "${SCALASCA_COMMAND} -instrument -comp=none -user ${CMAKE_CXX_COMPILER}"))
#        set(CMAKE_CXX_COMPILER "${SCALASCA_COMMAND} -instrument -comp=none -user ${CMAKE_CXX_COMPILER}")
#        endif(NOT ("${CMAKE_CXX_COMPILER}" STREQUAL "${SCALASCA_COMMAND} -instrument -comp=none -user ${CMAKE_CXX_COMPILER}"))
#        set(CXX "${SCALASCA_COMMAND} -instrument -comp=none -user ${CMAKE_CXX_COMPILER} ")

         # c) does not work, since RULE_LAUNCH_LINK also prefixes Archiver (ar)
#        set_property(GLOBAL PROPERTY RULE_LAUNCH_COMPILE "scalasca -instrument -comp=none -user ")
#        set_property(GLOBAL PROPERTY RULE_LAUNCH_LINK "scalasca -instrument -comp=none -user ")

        # add compile flags
		#add_cxx_flag("${SCALASCA_USER_CFLAGS}") -- scalasca wrapper does it itself
    	add_definitions(-DUG_PROFILER_SCALASCA)    
    
    # Vampir
    elseif("${PROFILER}" STREQUAL "Vampir")
        find_package(VampirTrace)
        if(VAMPIRTRACE_FOUND)
            message("-- Info: Vampir: using compiler wrapper c++: ${VAMPIRTRACE_CXX}")
            message("-- Info: Vampir: using compiler wrapper cc : ${VAMPIRTRACE_CC}")
            message("-- Info: Vampir: using inlcude dir: ${VAMPIRTRACE_INCLUDE_DIR}")
            message("-- Info: Vampir: using library dir: ${VAMPIRTRACE_LIBRARIES}")
            message("-- Info: Vampir: check that compiler wrapper are set!")
            message("--       Use: CC=\"vtcc\" CXX=\"vtcxx\" cmake ...")
            message("--       If not used: remove build completely and rerun cmake with prefix.")
        else(VAMPIRTRACE_FOUND)
        	message(FATAL_ERROR "PROFILER: ${PROFILER}: Cannot find required "
        	        "binary vtcxx and/or library or include paths. Make sure "
        	        "that PATH contains vtcxx executable.")
        endif(VAMPIRTRACE_FOUND)
    	add_definitions(-DUG_PROFILER_VAMPIR)    
    	add_definitions(-vt:inst manual -DVTRACE)

    # ScoreP
    elseif("${PROFILER}" STREQUAL "ScoreP")
        find_program(SCOREP_COMMAND scorep)
        if(SCOREP_COMMAND)
            message("-- Info: ScoreP: using compiler wrapper: ${SCOREP_COMMAND}")
            message("-- Info: ScoreP: check that compiler wrapper are set!")
            message("--       Use: CC=\"scorep --user mpicc\" CXX=\"scorep --user mpicxx\" cmake ...")
            message("--       If not used: remove build completely and rerun cmake with prefix.")
        else(SCOREP_COMMAND)
        	message(FATAL_ERROR "PROFILER: ${PROFILER}: Cannot find required "
        	        "binary scorep and/or library or include paths. Make sure "
        	        "that PATH contains scorep executable.")
        endif(SCOREP_COMMAND)
        
    	add_definitions(-DUG_PROFILER_SCOREP)    
    	#add_definitions(-vt:inst manual -DVTRACE)

    # wrong string in compiler
    else("${PROFILER}" STREQUAL "Shiny")
    	message(FATAL_ERROR "Unsupported PROFILER: ${PROFILER}. Options are: ${profilerOptions}")
    endif("${PROFILER}" STREQUAL "Shiny")
    
    # if one profiler is used this flag will be set
	add_definitions(-DUG_PROFILER)    # add to c++ flags
	set(UG_PROFILER ON)               # add Cmake variable
	
endif(NOT "${PROFILER}" STREQUAL "None")

########################################
# PROFILE_PCL
if(PROFILE_PCL)
	add_definitions(-DPROFILE_PCL)
endif(PROFILE_PCL)

# PROFILE_BRIDGE
if(PROFILE_BRIDGE)
	add_definitions(-DPROFILE_BRIDGE)
endif(PROFILE_BRIDGE)

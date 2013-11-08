# included from ug_includes.cmake
########################################
# PROFILER

if((NOT "${PROFILER}" STREQUAL "Shiny") AND SHINY_CALL_LOGGING)
    message(FATAL_ERROR " Shiny Call Logging activated but not Shiny. Use cmake -DPROFILER=Shiny ..")
endif((NOT "${PROFILER}" STREQUAL "Shiny") AND SHINY_CALL_LOGGING)


if(NOT("${PROFILER}" STREQUAL "None"))
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
	
endif(NOT("${PROFILER}" STREQUAL "None"))

########################################
# PROFILE_PCL
if(PROFILE_PCL)
	add_definitions(-DPROFILE_PCL)
endif(PROFILE_PCL)

# PROFILE_BRIDGE
if(PROFILE_BRIDGE)
	add_definitions(-DPROFILE_BRIDGE)
endif(PROFILE_BRIDGE)

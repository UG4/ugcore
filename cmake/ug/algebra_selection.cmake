# included from ug_includes.cmake
########################################
# CPU
# The cpu option sets defines for C and C++
if(GPU_ALGEBRA)
	MESSAGE(STATUS "Info: Using GPU Algebra.")
    add_definitions(-DUG_GPU)
endif()
if(CPU_ALGEBRA)
    MESSAGE(STATUS "Info: Using CPU Algebra.")
    if("${CPU}" STREQUAL "ALL")
        # todo checks for 4, VAR (-DUG_CPU_4 -DUG_CPU_VAR)!
        add_definitions(-DUG_CPU_1 -DUG_CPU_2 -DUG_CPU_3)

    else("${CPU}" STREQUAL "ALL")
        # CPU is a string of numbers (e.g. "1;2")
        # loop dims
        foreach(d ${CPU})
            # check if dim is valid
            if(d GREATER 4 OR d LESS 1)
                message(FATAL_ERROR "ERROR: Cannot build cpu blocksize ${d}. "
                                    "Valid options are: ${cpuOptions}")
            endif(d GREATER 4 OR d LESS 1)

            add_definitions(-DUG_CPU_${d})
        endforeach(d)
    endif("${CPU}" STREQUAL "ALL")
endif()
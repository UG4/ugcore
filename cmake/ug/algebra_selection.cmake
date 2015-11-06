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
# CPU
# The cpu option sets defines for C and C++
if(GPU_ALGEBRA)
	MESSAGE(STATUS "Info: Using GPU Algebra.")
    add_definitions(-DUG_GPU)
endif()
if(CPU_ALGEBRA)
    MESSAGE(STATUS "Info: Using CPU Algebra.")
    if("${CPU}" STREQUAL "ALL")
        # todo checks for 4, VAR (-DUG_CPU_4 -DUG_CPU_5 -DUG_CPU_VAR)!
        # todo: This is somehow misleading "ALL" != all posibilities, but only
        #       world 1-3. Should we fix this?!
        add_definitions(-DUG_CPU_1 -DUG_CPU_2 -DUG_CPU_3)

    else("${CPU}" STREQUAL "ALL")
        # CPU is a string of numbers (e.g. "1;2")
        # loop dims
        foreach(d ${CPU})
            # check if dim is valid
            if(d GREATER 5 OR d LESS 1)
                message(FATAL_ERROR "ERROR: Cannot build cpu blocksize ${d}. "
                                    "Valid options are: ${cpuOptions}")
            endif(d GREATER 5 OR d LESS 1)

            add_definitions(-DUG_CPU_${d})
        endforeach(d)
    endif("${CPU}" STREQUAL "ALL")
endif()
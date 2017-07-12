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
if(USE_LUAJIT)

 	#
	# TODO: Adjust paths here!
	# (or write automatic adjustment)
	# 
	SET(LUAJIT_INCLUDE_PATH "/opt/local/include/include/luajit-2.0")
	SET(LUAJIT_LIBRARY_PATH "/opt/local/lib/")
	SET(LUAJIT_LIBRARIES "luajit-5.1")
	
    if(STATIC_BUILD)
    	MESSAGE(STATUS "Info: LUAJIT requested, but static build. LUAJIT disabled.")
    	set(USE_LUAJIT OFF)
    else(STATIC_BUILD)
    
    	MESSAGE(STATUS "Info: Using LUAJIT, (include: ${LUAJIT_INCLUDE_PATH}, lib: ${LUAJIT_LIBRARY_PATH})")
    	
    	add_definitions(-DUSE_LUAJIT)
    	include_directories(${LUAJIT_INCLUDE_PATH})
    	link_directories(${LUAJIT_LIBRARY_PATH})
    	set(linkLibraries ${linkLibraries} ${LUAJIT_LIBRARIES})
		
		 # Make sure to use correct address space for luajit on Mac OS X
		 # (cf. 'Embedding LuaJIT' in luajit-2.0/doc/install.html)
		IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    		set(CMAKE_CXX_LINK_EXECUTABLE "${CMAKE_CXX_LINK_EXECUTABLE} -pagezero_size 10000 -image_base 100000000")
		ENDIF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
		
    endif(STATIC_BUILD)
endif(USE_LUAJIT)
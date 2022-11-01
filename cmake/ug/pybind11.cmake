# Copyright (c) 2020:  G-CSC, Goethe University Frankfurt
# Author: Arne Naegel
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


if(USE_PYBIND11)
    if(STATIC_BUILD)
    	MESSAGE(STATUS "Info: Pybind11 requested, but static build. Pybind11 disabled?")
    	SET(USE_PYBIND11 OFF)
    else(STATIC_BUILD)
    	MESSAGE(STATUS "Info: Using Pybind11")
    	
    	# Automatic
    	SET(pybind11_DIR ${UG_ROOT_CMAKE_PATH}/../../externals/pybind11/include)
    	
    	IF (EXTERNAL_PYBIND11)
    		# Automatic
    		FIND_PACKAGE(pybind11 CONFIG REQUIRED)
   	   	ENDIF(EXTERNAL_PYBIND11)
   	   	
   	   	MESSAGE(STATUS "Info: Using JSON from ${my_nlohmann_json_incl}")
    	include_directories(${my_nlohmann_json_incl})

    	add_definitions(-DUG_USE_PYBIND11)
    endif(STATIC_BUILD)
else(USE_PYBIND11)
	set(USE_PYBIND11 OFF)
endif(USE_PYBIND11)

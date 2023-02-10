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
if(USE_JSON)
    if(STATIC_BUILD)
    	MESSAGE(STATUS "Info: JSON requested, but static build. JSON disabled?")
    	SET(USE_JSON OFF)
    else(STATIC_BUILD)
    	MESSAGE(STATUS "Info: Using JSON")
    	
    	# Automatic
    	FIND_PACKAGE(nlohmann_json QUIET)
    	MESSAGE("-- Adding JSON from ${UG_ROOT_CMAKE_PATH}/../../externals/JSONForUG4/json-cxx/include")
    	include_directories(${UG_ROOT_CMAKE_PATH}/../../externals/JSONForUG4/json-cxx/include)
    	MESSAGE("-- Dir: ${NLOHMANN_JSON_INCLUDE_INSTALL_DIR}") 
    	add_definitions(-DUG_JSON)
    endif(STATIC_BUILD)
else(USE_JSON)
	set(USE_JSON OFF)
endif(USE_JSON)
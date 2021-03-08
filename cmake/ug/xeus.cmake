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


if(USE_XEUS) # included from ug_includes.cmake
    if(STATIC_BUILD)
    	MESSAGE(STATUS "Info: XEUS requested, but static build. XEUS disabled?")
    	SET(USE_XEUS OFF)
    else(STATIC_BUILD)
    	MESSAGE(STATUS "Info: Using XEUS")
    	
    	
		# Various Xeus libs
		find_package(xwidgets REQUIRED) 
		find_package(xeus REQUIRED)
		find_package(xtl REQUIRED)

		message(STATUS "Jupyter-Plugin: xwidgets_FOUND=${xwidgets_FOUND} ${xwidgets_INCLUDE_DIRS}")
		message(STATUS "Jupyter-Plugin: xeus_FOUND=${xeus_FOUND} ${xeus_INCLUDE_DIRS} ")
		message(STATUS "Jupyter-Plugin: xtl_FOUND=${xtl_FOUND} ${xtl_INCLUDE_DIRS}")


		# Find libraries
		if(NOT DEFINED ${xeus_LIBRARY})
			find_library(xeus_LIBRARY NAMES xeus PATHS "${xtl_INCLUDE_DIRS}/../lib")
		endif() 

		if(NOT DEFINED ${xwidgets_LIBRARY})
			find_library(xwidgets_LIBRARY NAMES xtl PATHS "${xtl_INCLUDE_DIRS}/../lib")
		endif() 

		message(STATUS "Jupyter-Plugin: using  ${xeus_LIBRARY} ${xwidgets_LIBRARY}")   




		# We also need  'zeromq' and 'sodium'
		## load in pkg-config support
		set(ENV{PKG_CONFIG_PATH} "$ENV{CONDA_PREFIX}/lib/pkgconfig")
		message(STATUS "Jupyter-Plugin: Checking pkg-config $ENV{PKG_CONFIG_PATH}")
		
		find_package(PkgConfig REQUIRED)
		pkg_search_module(ZEROMQ REQUIRED libzmq cppzmq zmq) ## use pkg-config to get hints for 0mq locations
		pkg_search_module(SODIUM REQUIRED libsodium sodium)  ## use pkg-config to get hints for sodium locations

		message(STATUS "Jupyter-Plugin:${ZEROMQ_FOUND} ${ZEROMQ_LINK_LIBRARIES} ${ZEROMQ_LIBRARIES} ")
		message(STATUS "Jupyter-Plugin:${SODIUM_FOUND} ${SODIUM_LINK_LIBRARIES} ${SODIUM_LIBRARIES}")
		
		## use the hint from about to find the location of libzmq, sodium
		if(NOT DEFINED ${ZEROMQ_LIBRARY_DIRS})
			message(WARNING "Jupyter-Plugin: $ZEROMQ_LIBRARY_DIRS not defined.")
			set (PC_ZeroMQ_LIBRARY_DIRS "${xtl_INCLUDE_DIRS}/../lib")
		endif()	

		# find_library(ZeroMQ_LIBRARY NAMES zmq PATHS ${PC_ZeroMQ_LIBRARY_DIRS})
		# find_library(sodium_LIBRARY NAMES sodium PATHS ${PC_ZeroMQ_LIBRARY_DIRS})
    	
    	
    	# Automatic
    	# FIND_PACKAGE(nlohmann_json QUIET)
    	# MESSAGE("-- Adding JSON from ${UG_ROOT_CMAKE_PATH}/../../externals/JSONForUG4/json-cxx")	
    	# include_directories(${UG_ROOT_CMAKE_PATH}/../../externals/JSONForUG4/json-cxx)
    	# MESSAGE("-- Dir: ${NLOHMANN_JSON_INCLUDE_INSTALL_DIR}") 
    	add_definitions(-DUG_XEUS)
    endif(STATIC_BUILD)
else(USE_XEUS)
	set(USE_XEUS OFF)
endif(USE_XEUS)
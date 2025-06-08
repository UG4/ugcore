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

message(STATUS "This is pybind11.cmake")
if(USE_PYBIND11)

cmake_minimum_required(VERSION 3.12)

# ug4pybind_add_module
# argument1: name of Python module
# argument2: source files
# argument3: source files
message(STATUS "Installing ug4pybind_add_module")
function(ug4pybind_add_module pyPluginName myPluginSources myLinkLibs)
  	
  	# Add module.
	pybind11_add_module(${pyPluginName} ${myPluginSources})
	set_target_properties(${pyPluginName} PROPERTIES	
  		CXX_STANDARD 11 
  		CXX_STANDARD_REQUIRED YES
   		CXX_EXTENSIONS NO
   		OUTPUT_NAME ug4py/${pyPluginName})
   		
   	# Link against 'this' plugin and UG4 library.
	target_link_libraries (${pyPluginName} PRIVATE ${myLinkLibs})
    
endfunction(ug4pybind_add_module)


# get_libug4_from_pip
# return_value1: location of library in current pip directory.
function(get_libug4_from_pip RETURN_VAR)
	SET(UG4_STATIC_LIB_BASENAME ug4_s)
	#if(WIN32)
    # 	set(UG4_STATIC_LIB_NAME "${LIB_BASENAME}.lib")
	#else()
    #	set(UG4_STATIC_LIB_NAME "lib${LIB_BASENAME}.a")
	#endif()
	
	# Extract possible location of pip-package
	execute_process(
    		#COMMAND pip show ug4py-base
			COMMAND python3 -m pip show ug4py-base
    		OUTPUT_VARIABLE PIP_SHOW_OUTPUT
    		OUTPUT_STRIP_TRAILING_WHITESPACE
	)

	# Extract the "Location" line & remove the ' prefix to get the actual path
	string(REGEX MATCH "Location: .+" LOCATION_LINE "${PIP_SHOW_OUTPUT}")
	string(REPLACE "Location: " "" UG4PY_BASE_PACKAGE_PATH "${LOCATION_LINE}")
	string(REGEX MATCH "^[^\n\r]*" UG4PY_BASE_PACKAGE_PATH "${UG4PY_BASE_PACKAGE_PATH}") # remove lines.
	message(STATUS "Package path: '${UG4PY_BASE_PACKAGE_PATH}'")
	
	# Find library.
	find_library(PYUG4_STATIC_LIBRARY
    			NAMES ${UG4_STATIC_LIB_BASENAME} lib${UG4_STATIC_LIB_BASENAME}
    			PATHS ${UG4PY_BASE_PACKAGE_PATH}/lib ${UG4PY_BASE_PACKAGE_PATH}/lib64)
	
	message(STATUS "Library ${UG4_STATIC_LIB_BASENAME} found: ${PYUG4_STATIC_LIBRARY}")
	
	set(${RETURN_VAR} "${PYUG4_STATIC_LIBRARY}" PARENT_SCOPE) 
endfunction()

endif(USE_PYBIND11)


if(USE_PYBIND11)
	# Check for Python.
	MESSAGE(STATUS "Info: ****************************************************************************************")
    MESSAGE(STATUS "Info: Pybind11 enabled.")
 	
 	# Check for Python.
    # FIND_PACKAGE(Python REQUIRED COMPONENTS Interpreter Development)
	FIND_PACKAGE(Python REQUIRED COMPONENTS Interpreter Development.Module)
    MESSAGE(STATUS "Info: Found Python = ${Python_FOUND}")
    MESSAGE(STATUS "Info: Using INC${Python_INCLUDE_DIRS}")
    MESSAGE(STATUS "Info: Using LIB ${Python_LIBRARIES}")
	MESSAGE(STATUS "Info: Using ID=${Python_INTERPRETER_ID}")
	MESSAGE(STATUS "Info: Using VERSION=${Python_VERSION}")
    	
    # Then check for Pybind11 (order is important!).
	add_subdirectory(
		${UG_ROOT_CMAKE_PATH}/../../externals/PybindForUG4/pybind11 # Source dir.
		{UG_ROOT_CMAKE_PATH}/bin   									# Binary dir (mandatory, but not used).
	)
   	MESSAGE(STATUS "Info: Found Pybind11 = ${pybind11_FOUND}")
   	MESSAGE(STATUS "Info: Using Pybind11 INC ${pybind11_INCLUDE_DIR}")
   	MESSAGE(STATUS "Info: Using Pybind11 INC ${pybind11_INCLUDE_DIRS}")
   	MESSAGE(STATUS "Info: ****************************************************************************************")
    	
	# Change policy in order to avoid errors.
    cmake_policy(SET CMP0057 NEW)

	# Expand includes.
    include_directories(${pybind11_INCLUDE_DIRS})
    	
	# Expand libraries.
	set(linkLibraries ${linkLibraries} ${Python_LIBRARIES})
	
	# Set define for UG4.
	add_definitions(-DUG_USE_PYBIND11)

	# Extend flags
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fPIC -fvisibility=hidden")

    if(STATIC_BUILD)
		# For Pybind with static builds, embedded plugins are disabled.
		# Instead, individual static libs are build as temporaries
		# and the Python-API lib are linked against these.
    	MESSAGE(STATUS "Info: Pybind11 w/ static build. Embedded plugins are disabled.")
		
		# This undos settings in ug_include.cmake
		set(buildEmbeddedPlugins OFF)
    endif(STATIC_BUILD)

else(USE_PYBIND11)
	set(USE_PYBIND11 OFF)
endif(USE_PYBIND11)

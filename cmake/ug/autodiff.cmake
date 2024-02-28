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



# Set Autodiff include path
if(USE_AUTODIFF)
  	# This sections adds autodiff to the include path.
	# Note: Autodiff requires C++17-compliant compilers

	IF (EXTERNAL_AUTODIFF) # Automatic
    	FIND_PACKAGE(autodiff CONFIG REQUIRED)
	else (EXTERNAL_AUTODIFF)# Builtin
		SET(autodiff_DIR ${UG_ROOT_CMAKE_PATH}/../../externals/AutodiffForUG4/autodiff)
   	ENDIF(EXTERNAL_AUTODIFF)
   	   	
   	MESSAGE(STATUS "Info: Using Autodiff from ${autodiff_DIR}")
	MESSAGE(STATUS "Info: ... compiling plugins will require at least C++17.")
     
	# These are the important lines...
	include_directories(${autodiff_DIR})
    add_definitions(-DUG_USE_AUTODIFF)


	# As a next step, update the submodule.
	find_package(Git QUIET)
	if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
		# Update submodules as needed
		option(GIT_SUBMODULE "Check submodules during build" ON)
		if(GIT_SUBMODULE)
			message(STATUS "Submodule update:${GIT_EXECUTABLE}  ${autodiff_DIR}/..")
			execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
							WORKING_DIRECTORY "${autodiff_DIR}/.."
							RESULT_VARIABLE GIT_SUBMOD_RESULT)
			if(NOT GIT_SUBMOD_RESULT EQUAL "0")
				message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
			endif()
		endif()
	endif()
	
	if(NOT EXISTS "${autodiff_DIR}/CMakeLists.txt")
		message(ERROR "The submodules in ${autodiff_DIR} were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
	endif()


else(USE_AUTODIFF)
	set(USE_AUTODIFF OFF)
endif(USE_AUTODIFF)



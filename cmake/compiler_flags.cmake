# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Martin Scherer
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

# This function sets the CMAKE_C_FLAGS and CMAKE_CXX_FLAGS variabels.
# 
# TODO: maybe add sanity checks for flags
#

# add flag for c language only
function (add_c_flag flag)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${flag}" CACHE STRING "overriden flags!" FORCE)
endfunction(add_c_flag)

# add flag for c++ language only
function (add_cpp_flag flag)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${flag}" CACHE STRING "overriden flags!" FORCE)
endfunction(add_cpp_flag)

# add flag for c and c++ language
function(add_cxx_flag flag)
	add_c_flag(${flag})
	add_cpp_flag(${flag})
endfunction(add_cxx_flag)

################################################################################
# Used to clear c/cxx flags. Handles environment variables CXX_FLAGS and C_FLAGS
function(reset_cxx_flags)
	# check env for set flags, and overwrite present cxx flag string
	# note: this allows also empty strings "" 
	if(DEFINED ENV{CXX_FLAGS})
		set(CXX_FLAGS "$ENV{CXX_FLAGS}" CACHE STRING "custom cxx flags" FORCE)
	endif()
	
	if(DEFINED ENV{C_FLAGS})
		set(C_FLAGS "$ENV{C_FLAGS}" CACHE STRING "custom c flags" FORCE)
	endif()
	
	# first reset flags, then check env for given flags
	foreach(lang C CXX)
		set(CMAKE_${lang}_FLAGS "" CACHE STRING "clear flags" FORCE)
	endforeach()
	
	# re-add custom flags for c and c++
	if(CXX_FLAGS)
		add_cpp_flag(${CXX_FLAGS})
	endif()
	
	if(C_FLAGS)
		add_c_flag(${C_FLAGS})
	endif()
endfunction(reset_cxx_flags)

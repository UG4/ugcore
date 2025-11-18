# Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter, Andreas Vogel
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

################################################
# ugbase
################################################

cmake_minimum_required(VERSION 2.8.12...3.27.1)

project(P_UGBASE)

# set global output paths for UGBASE
# They are relative to the source dir for which cmake was executed.
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

set(Sources	ug.cpp)

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase ${Sources})
else(BUILD_ONE_LIB)
	add_library(ugbase ${Sources})
endif(BUILD_ONE_LIB)

########################
# information on the current compilation environment
add_subdirectory(compile_info)

########################
# common
add_subdirectory(common)

########################
# lib_grid
if(buildGrid)
	add_subdirectory(lib_grid)
endif(buildGrid)

########################
# lib_registry
if(buildRegistry)
	add_subdirectory(registry)
endif(buildRegistry)

########################
# bindings
if(buildBindings)
	add_subdirectory(bindings/lua)
	if(buildForVRL)
    	add_subdirectory(bindings/vrl)
	endif(buildForVRL)
endif(buildBindings)

########################
# lib_algebra
if(buildAlgebra)
	add_subdirectory(lib_algebra)	
endif(buildAlgebra)

########################
# lib_disc
if(buildDisc)
	add_subdirectory(lib_disc)
endif(buildDisc)

########################
# pcl
if(PARALLEL)
	add_subdirectory(pcl)
endif(PARALLEL)

########################
# bridge
if(buildBridge)
	add_subdirectory(bridge)
endif(buildBridge)

########################
# ug4 shell
if(buildUGShell)
    add_subdirectory(ug_shell)
endif(buildUGShell)

if(INTERNAL_BOOST)
	add_subdirectory(../../externals/BoostForUG4/ externals/BoostForUG4/)
endif(INTERNAL_BOOST)
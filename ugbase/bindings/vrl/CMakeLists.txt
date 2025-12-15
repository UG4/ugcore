# Copyright (c) 2010-2014:  Steinbeis Forschungszentrum (STZ Ölbronn)
# Author: Michael Hoffer
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

cmake_minimum_required(VERSION 3.10)

####
# bindings_vrl Library
####

project(P_BINDINGS_VRL)

set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/lib)
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/bin)

include_directories(${P_UGBASE_SOURCE_DIR})

include("../../../cmake/ug_includes.cmake")


set(Sources	bindings_vrl.cpp
			type_converter.cpp
			canvas.cpp
			messaging.cpp
			invocation.cpp
			playground.cpp
			user_data.cpp
			basicTest.cpp
)

if(buildAlgebra)
	set(Sources ${Sources} vrl_bridge.cpp)
endif(buildAlgebra)


if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/bindings/vrl  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bindings_vrl  ${Sources})
endif(BUILD_ONE_LIB)

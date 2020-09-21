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
# C++11
IF(CXX11)
	IF(${CMAKE_CXX_COMPILER_ID} MATCHES GNU|Clang|Intel)
		# Check for the compilers's C++11 capabilities
		INCLUDE(CheckCXXCompilerFlag)
		CHECK_CXX_COMPILER_FLAG(-std=c++0x HAVE_CXX0X)
		# since GCC4.7 (c++0x will be removed in future versions of GCC)
		CHECK_CXX_COMPILER_FLAG(-std=c++11 HAVE_CXX11)
		# Add appropriate compiler flags
		IF(HAVE_CXX11)
			SET(CXX11_FLAG "-std=c++11")
		ELSEIF(HAVE_CXX0X)
			SET(CXX11_FLAG "-std=c++0x")
		ENDIF()

		IF(CXX11_FLAG)
			ADD_DEFINITIONS(-DUG_CXX11)
			add_cpp_flag(${CXX11_FLAG})
			MESSAGE(STATUS "Info: C++11 enabled. (flag: ${CXX11_FLAG})")
		ELSE()
			SET(CXX11 OFF)
			MESSAGE(STATUS "Info: Compiler does not support C++11 standard.")
		ENDIF()
	ELSE()
		MESSAGE(STATUS "Info: Enabling C++11 is currently only supported with GCC, Clang, or Intel")
	ENDIF()
ENDIF(CXX11)

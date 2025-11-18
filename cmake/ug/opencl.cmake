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
if (OpenCL)
	if (APPLE)
		message (STATUS "Info: Found OpenCL: -framework OpenCL (Apple)")
		set (OPENCL_FOUND YES)
		set (linkLibraries "${linkLibraries} -framework OpenCL")
		add_definitions (-DOPENCL_AVAILABLE)

	else ()
		INCLUDE ("${UG_ROOT_PATH}/cmake/FindOpenCL.cmake")	
		if (OPENCL_FOUND)
			message (STATUS "Info: Found OpenCL. Inc: ${OPENCL_INCLUDE_DIRS}, lib: ${OPENCL_LIB_DIR}")
			include_directories (${OPENCL_INCLUDE_DIRS})
			set (linkLibraries ${linkLibraries} ${OPENCL_LIB_DIR})
			add_definitions (-DOPENCL_AVAILABLE)

		else ()
			message (FATAL_ERROR " OpenCL not found! Aborting.")
		endif ()
	endif ()
endif ()
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
if (CUDA)
	set (MY_HOSTNAME $ENV{HOSTNAME})
	if (MY_HOSTNAME STREQUAL "tales")
		set (CUDA_TOOLKIT_ROOT_DIR /usr/local/cuda-5.0)
    endif ()
    find_package (CUDA)
	if (CUDA_FOUND)
		set (CUDA_NVCC_FLAGS "-arch=sm_13")
		message (STATUS "Info: CUDA ${CUDA_VERSION}. Toolkit root dir: ${CUDA_TOOLKIT_ROOT_DIR}, nvcc: ${CUDA_NVCC_EXECUTABLE}")

		file (WRITE ${CMAKE_BINARY_DIR}/cuda_jit.h
			"#define CUDA_TOOLKIT_ROOT_DIR \"${CUDA_TOOLKIT_ROOT_DIR}\"\n"
			"#define CUDA_NVCC_EXECUTABLE \"${CUDA_NVCC_EXECUTABLE}\"\n"
			"#define CUDA_VERSION \"${CUDA_VERSION}\"\n"
			"#define CUDA_NVCC_FLAGS \"${CUDA_NVCC_FLAGS}\"\n")

		include_directories (${CUDA_TOOLKIT_INCLUDE})
		include_directories (${CUDA_TOOLKIT_ROOT_DIR}/samples/common/inc)

		set (linkLibraries ${linkLibraries} ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY})
		## cusparse is not found on some systems
		#set(linkLibraries ${linkLibraries} "${CUDA_TOOLKIT_ROOT_DIR}/lib64/libcusparse.so")
		set (linkLibraries ${linkLibraries} ${CUDA_cusparse_LIBRARY})
		add_definitions (-DCUDA_AVAILABLE)
	else (CUDA_FOUND)
		message (FATAL_ERROR "Error: Couldn't find CUDA")
	endif ()
endif ()
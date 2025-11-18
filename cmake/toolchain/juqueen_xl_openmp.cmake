# Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
# Author: Stephan Grein
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

################################################################################
# created by Stephan Grein <stephan.grein@gcsc.uni-frankfurt.de>
# derived from juqueen.cmake established by Ingo Heppner.
#
# enables OpenMP for IBM XL compiler and associated tools on JUQUEEN (BlueGene/Q)
################################################################################

# set CMAKE_SYSTEM_NAME to include automatically corresponding platform files
set (CMAKE_SYSTEM_NAME BlueGeneQ-static)
list (APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules")

# no dynamic build
set (enableDynamicOption OFF)

# add -fno-strict-aliasing option to the compiler flags
set (enableNoStrictAliasingOption ON)

# variables (note _r suffix indicates OpenMP compiler)
set (GCC_ROOT  "/bgsys/drivers/ppcfloor/gnu-linux")
set (GCC_NAME  "powerpc64-bgq-linux")
set (MPI_ROOT  "/bgsys/drivers/ppcfloor/comm/gcc")

# serial GCC
set (CMAKE_C_COMPILER       ${GCC_ROOT}/bin/${GCC_NAME}-gcc)
set (CMAKE_CXX_COMPILER     ${GCC_ROOT}/bin/${GCC_NAME}-g++)
set (CMAKE_Fortran_COMPILER ${GCC_ROOT}/bin/${GCC_NAME}-gfortran)

# mpi wrappers for GCC (note _r suffix indicates OpenMP compiler)
set (MPI_C_COMPILER       ${MPI_ROOT}/bin/mpicc)
set (MPI_CXX_COMPILER     ${MPI_ROOT}/bin/mpicxx)
set (MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpif90)

message (STATUS "TMP INFO: Value of '\${CMAKE_C_COMPILER}'       is: ${CMAKE_C_COMPILER}")       # TMP
message (STATUS "TMP INFO: Value of '\${CMAKE_CXX_COMPILER}'     is: ${CMAKE_CXX_COMPILER}")     # TMP
message (STATUS "TMP INFO: Value of '\${CMAKE_Fortran_COMPILER}' is: ${CMAKE_Fortran_COMPILER}") # TMP

message (STATUS "TMP INFO: Value of '\${MPI_C_COMPILER}'         is: ${MPI_C_COMPILER}")         # TMP
message (STATUS "TMP INFO: Value of '\${MPI_CXX_COMPILER}'       is: ${MPI_CXX_COMPILER}")       # TMP
message (STATUS "TMP INFO: Value of '\${MPI_Fortran_COMPILER}'   is: ${MPI_Fortran_COMPILER}")   # TMP

# add compiler flag -fopenmp (OpenMP)
set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -fopenmp")
set ( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -fopenmp")

# XL specific options (not known by GCC to my knowledge)
set ( CMAKE_C_FLAGS  "${CMAKE_C_FLAGS} -qsmp=omp -qnosave")
set ( CMAKE_CXX_FLAGS  "${CMAKE_CXX_FLAGS} -qsmp=omp -qnosave")

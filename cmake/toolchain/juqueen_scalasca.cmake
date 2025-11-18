# Copyright (c) 2013:  G-CSC, Goethe University Frankfurt
# Author: Ingo Heppner, Andreas Vogel
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
# created by Ingo Heppner
# Ingo.Heppner@gcsc.uni-frankfurt.de
#
# Toolchain file for JuQueen (type IBM Blue Gene/Q) FZ Juelich.
#
# This toolchain file is included in a CMake run by (check with 'cmake --trace')
# 'Modules/CMakeDetermineSystem.cmake', and (later again) by the - generated -
# file '<ug4_build_dir>/CMakeFiles/CMakeSystem.cmake'.
# (Subdirectory 'Modules/' is part of your CMake installation, on JuQueens login
# nodes: '/usr/local/cmake/share/cmake-2.8/Modules/').
#
# Appropriate platform files for the "CMake system name" and compilers chosen
# here are called later (by 'Modules/CMakeSystemSpecificInformation.cmake'):
#
# I.e., for "BlueGeneQ-static" as "CMake system name" the platform file
# 'Modules/Platform/BlueGeneQ-static.cmake' is included first, which includes
# 'Modules/Platform/BlueGeneQ-base.cmake'. Afterwards CMake calls
# 'Modules/Platform/BlueGeneQ-static-GNU-C.cmake' and
# 'Modules/Platform/BlueGeneQ-static-GNU-CXX.cmake'
# if GNU compilers are chosen (by 'Modules/CMakeCInformation.cmake' and
# 'Modules/CMakeCXXInformation.cmake' respectively).
#
################################################################################
#
# Attention:
# For this to work a CMake version with platform files for Blue Gene/Q is required!
# The currently (October 2012) installed version on JuQueen does not (the FZJ
# support is asked to support one).
#
# In the meantime one can resort to a local installation of (the most recent
# version of) CMake, where appropriate platform files are copied into the
# 'Modules/Platform/' sub directory:
#
#  ~/local/share/cmake-2.8/Modules/Platform/BlueGeneQ-*.cmake
#
# The following is based on information from
# 'https://code.google.com/p/elemental/source/browse/cmake/toolchains/BGQ-gnu-netlib.cmake',
# utilising platform files from
# 'http://www.cmake.org/Bug/bug_relationship_graph.php?bug_id=13512&graph=dependency'
# (found 03092012), which are copied into the 'Modules/Platform/' directory of
# a locally installed CMake, version 2.8.9 (installed 03092012).
#
################################################################################

# Important: Setting the "cmake system name" will lead to automatic inclusion of
# the corresponding platform files:
set (CMAKE_SYSTEM_NAME BlueGeneQ-static)

# This option tells cmake/ug_includes.cmake to add the -dynamic option to the compiler flags.
set (enableDynamicOption OFF)

# This option tells cmake/ug_includes.cmake to add the -fno-strict-aliasing option to the compiler flags.
set (enableNoStrictAliasingOption ON)

# Some variables for pathes
set (GCC_ROOT  "/bgsys/drivers/ppcfloor/gnu-linux")
set (GCC_NAME  "powerpc64-bgq-linux")
set (MPI_ROOT  "/bgsys/drivers/ppcfloor/comm/gcc")
#set (PAMI_ROOT "/bgsys/drivers/ppcfloor/comm/sys")
#set (SPI_ROOT  "/bgsys/drivers/ppcfloor/spi")

# The serial GNU compilers
set (CMAKE_C_COMPILER       "scalasca -instrument -comp=none -user ${GCC_ROOT}/bin/${GCC_NAME}-gcc")
set (CMAKE_CXX_COMPILER     "scalasca -instrument -comp=none -user ${GCC_ROOT}/bin/${GCC_NAME}-g++")
set (CMAKE_Fortran_COMPILER "scalasca -instrument -comp=none -user ${GCC_ROOT}/bin/${GCC_NAME}-gfortran")

# The MPI wrappers for the GNU compilers
set (MPI_C_COMPILER       "scalasca -instrument -comp=none -user ${MPI_ROOT}/bin/mpicc")
set (MPI_CXX_COMPILER     "scalasca -instrument -comp=none -user ${MPI_ROOT}/bin/mpicxx")
set (MPI_Fortran_COMPILER "scalasca -instrument -comp=none -user ${MPI_ROOT}/bin/mpif90")

message (STATUS "TMP INFO: Value of '\${CMAKE_C_COMPILER}'       is: ${CMAKE_C_COMPILER}")       # TMP
message (STATUS "TMP INFO: Value of '\${CMAKE_CXX_COMPILER}'     is: ${CMAKE_CXX_COMPILER}")     # TMP
message (STATUS "TMP INFO: Value of '\${CMAKE_Fortran_COMPILER}' is: ${CMAKE_Fortran_COMPILER}") # TMP

message (STATUS "TMP INFO: Value of '\${MPI_C_COMPILER}'         is: ${MPI_C_COMPILER}")         # TMP
message (STATUS "TMP INFO: Value of '\${MPI_CXX_COMPILER}'       is: ${MPI_CXX_COMPILER}")       # TMP
message (STATUS "TMP INFO: Value of '\${MPI_Fortran_COMPILER}'   is: ${MPI_Fortran_COMPILER}")   # TMP

# For debugging purposes
#include(CMakePrintSystemInformation)

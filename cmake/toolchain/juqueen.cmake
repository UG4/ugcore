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
set(CMAKE_SYSTEM_NAME BlueGeneQ-static)
set(CMAKE_MODULE_PATH ${CMAKE_SOURCE_DIR}/cmake/platform}

# This option tells cmake/ug_includes.cmake to add the -dynamic option to the compiler flags.
SET(enableDynamicOption OFF)

# This option tells cmake/ug_includes.cmake to add the -fno-strict-aliasing option to the compiler flags.
SET(enableNoStrictAliasingOption ON)

# Some variables for pathes
set(GCC_ROOT  "/bgsys/drivers/ppcfloor/gnu-linux")
set(GCC_NAME  "powerpc64-bgq-linux")
set(MPI_ROOT  "/bgsys/drivers/ppcfloor/comm/gcc")
#set(PAMI_ROOT "/bgsys/drivers/ppcfloor/comm/sys")
#set(SPI_ROOT  "/bgsys/drivers/ppcfloor/spi")

# The serial GNU compilers
set(CMAKE_C_COMPILER       ${GCC_ROOT}/bin/${GCC_NAME}-gcc)
set(CMAKE_CXX_COMPILER     ${GCC_ROOT}/bin/${GCC_NAME}-g++)
set(CMAKE_Fortran_COMPILER ${GCC_ROOT}/bin/${GCC_NAME}-gfortran)

# The MPI wrappers for the GNU compilers
set(MPI_C_COMPILER       ${MPI_ROOT}/bin/mpicc)
set(MPI_CXX_COMPILER     ${MPI_ROOT}/bin/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpif90)

message(STATUS "TMP INFO: Value of '\${CMAKE_C_COMPILER}'       is: ${CMAKE_C_COMPILER}")       # TMP
message(STATUS "TMP INFO: Value of '\${CMAKE_CXX_COMPILER}'     is: ${CMAKE_CXX_COMPILER}")     # TMP
message(STATUS "TMP INFO: Value of '\${CMAKE_Fortran_COMPILER}' is: ${CMAKE_Fortran_COMPILER}") # TMP

message(STATUS "TMP INFO: Value of '\${MPI_C_COMPILER}'         is: ${MPI_C_COMPILER}")         # TMP
message(STATUS "TMP INFO: Value of '\${MPI_CXX_COMPILER}'       is: ${MPI_CXX_COMPILER}")       # TMP
message(STATUS "TMP INFO: Value of '\${MPI_Fortran_COMPILER}'   is: ${MPI_Fortran_COMPILER}")   # TMP

# For debugging purposes
#include(CMakePrintSystemInformation)

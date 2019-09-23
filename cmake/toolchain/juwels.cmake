################################################################################
# created by Markus Breit
# markus.breit@gcsc.uni-frankfurt.de
# adapted from corresponding JuQueen file created by Ingo Heppner
#
# Toolchain file for JUWELS (type Dual Intel Xeon Skylake 8168) FZ Juelich.
#
# This toolchain file is included in a CMake run by (check with 'cmake --trace')
# 'Modules/CMakeDetermineSystem.cmake', and (later again) by the - generated -
# file '<ug4_build_dir>/CMakeFiles/CMakeSystem.cmake'.
# (Subdirectory 'Modules/' is part of your CMake installation, on JURECA's login
# nodes: '/usr/local/software/jureca/Stage3/software/Core/CMake/3.2.3/share/cmake-3.2/Modules').
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

# Important: Setting the "cmake system name" will lead to automatic inclusion of
# the corresponding platform files:
set(CMAKE_SYSTEM_NAME Linux)
#list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/ugcore/cmake/modules")

# This option tells cmake/ug_includes.cmake to add the -dynamic option to the compiler flags.
SET(enableDynamicOption OFF)

# This option tells cmake/ug_includes.cmake to add the -fno-strict-aliasing option to the compiler flags.
SET(enableNoStrictAliasingOption ON)

# Some variables for paths
set(GCC_ROOT  "/gpfs/software/juwels/stages/2018b/software/icc/2019.0.117-GCC-7.3.0/compilers_and_libraries_2019.0.117/linux/")
set(MPI_ROOT  "/gpfs/software/juwels/stages/2018b/software/psmpi/5.2.1-1-iccifort-2019.0.117-GCC-7.3.0/")
set(FTR_ROOT  "/gpfs/software/juwels/stages/2018b/software/ifort/2019.0.117-GCC-7.3.0/compilers_and_libraries_2019.0.117/linux/")

# The serial GNU compilers
set(CMAKE_C_COMPILER       ${GCC_ROOT}/bin/intel64/icc)
set(CMAKE_CXX_COMPILER     ${GCC_ROOT}/bin/intel64/icpc)
set(CMAKE_Fortran_COMPILER ${FTR_ROOT}/bin/intel64/ifort)

# The MPI wrappers for the GNU compilers
set(MPI_C_COMPILER       ${MPI_ROOT}/bin/mpicc)
set(MPI_CXX_COMPILER     ${MPI_ROOT}/bin/mpicxx)
set(MPI_Fortran_COMPILER ${MPI_ROOT}/bin/mpif90)

#message(STATUS "TMP INFO: Value of '\${CMAKE_C_COMPILER}'       is: ${CMAKE_C_COMPILER}")       # TMP
#message(STATUS "TMP INFO: Value of '\${CMAKE_CXX_COMPILER}'     is: ${CMAKE_CXX_COMPILER}")     # TMP
#message(STATUS "TMP INFO: Value of '\${CMAKE_Fortran_COMPILER}' is: ${CMAKE_Fortran_COMPILER}") # TMP

message(STATUS "TMP INFO: Value of '\${MPI_C_COMPILER}'         is: ${MPI_C_COMPILER}")         # TMP
message(STATUS "TMP INFO: Value of '\${MPI_CXX_COMPILER}'       is: ${MPI_CXX_COMPILER}")       # TMP
message(STATUS "TMP INFO: Value of '\${MPI_Fortran_COMPILER}'   is: ${MPI_Fortran_COMPILER}")   # TMP

# For debugging purposes
#include(CMakePrintSystemInformation)

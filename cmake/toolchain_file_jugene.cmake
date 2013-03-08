# Based on the scripts from Ingo Heppner.

# This option tells ug_cmake_includes.txt to add the -dynamic option to the compiler flags.
SET(enableDynamicOption ON)

SET(CMAKE_C_COMPILER       /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc)
SET(CMAKE_CXX_COMPILER     /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-g++)
SET(CMAKE_Fortran_COMPILER /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gfortran)

include(Platform/BlueGeneP-base)

__BlueGeneP_set_dynamic_flags(GNU CXX)

set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")


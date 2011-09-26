# Based on the scripts from Ingo Heppner.

#SET(CMAKE_SYSTEM_NAME BlueGeneP)
#SET(CMAKE_SYSTEM_VERSION 1)

SET(CMAKE_C_COMPILER       /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gcc)
SET(CMAKE_CXX_COMPILER     /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-g++)
SET(CMAKE_Fortran_COMPILER /bgsys/drivers/ppcfloor/gnu-linux/bin/powerpc-bgp-linux-gfortran)

include(Platform/BlueGeneP-base)

__BlueGeneP_set_dynamic_flags(GNU CXX)

set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")

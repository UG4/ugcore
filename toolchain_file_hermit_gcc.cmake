# created by Martin Rupp
# martin.rupp@gcsc.uni-frankfurt.de

# toolchain file for Hermit/XE6 HLRS Stuttgart with GNU C Compiler

# USAGE:
# first, change programming environment:
#   module swap $(module li 2>&1 | awk '/PrgEnv/{print $2}') PrgEnv-gnu
# then execute cmake
#   cmake -DCMAKE_TOOLCHAIN_FILE=../toolchain_file_hermit_gcc.cmake ..

# check modules:
# module list
# also: module avail / module load / module swap.

# for use on the CRAY XE6 (Hermit) of HLRS/Karlsruhe
# http://www.hlrs.de/systems/platforms/cray-xe6-hermit/

# see also $CRAY_* environment variables
# $CRAY_UGNI_POST_LINK_OPTS

# this seems to be not necessary for GCC
# we leave it in since it could be that -rdynamic is slower
SET(CMAKE_SYSTEM_NAME Catamount)

SET(CMAKE_Fortran_COMPILER ftn)

SET(CMAKE_C_COMPILER cc CACHE FORCE "")
SET(CMAKE_CXX_COMPILER CC CACHE FORCE "")
SET(BUILTIN_LAPACK YES CACHE FORCE "")
SET(BUILTIN_BLAS YES CACHE FORCE "")
SET(BUILTIN_MPI YES CACHE FORCE "")

SET(STATIC ON CACHE FORCE "")


# Khabi stuff
#SET(MPI_INCLUDE_PATH $(MPICH_DIR)/include)
#SET(MPI_LIBRARY -L$(MPICH_DIR)/lib)
#set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
#set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
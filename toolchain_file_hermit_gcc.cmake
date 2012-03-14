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

# to use cmake 2.8.7:
# module load tools/cmake/2.8.7
IF(${CMAKE_VERSION} STRLESS "2.8.7")
	# emulating Hermit.cmake
	SET(CMAKE_SYSTEM_NAME Catamount)
	SET(MPI_Fortran_NO_INTERROGATE CMAKE_Fortran_COMPILER)
	SET(MPI_LIBRARY -L$(MPICH_DIR)/lib)
	SET(MPI_EXTRA_LIBRARY -L$(MPICH_DIR)/lib)
	SET(BLA_STATIC ON)
	SET(BLA_VENDOR All)
	SET(BLAS_FIND_QUIETLY ON)
	SET(LAPACK_FIND_QUIETLY ON)
	SET(LAPACK_LIBRARIES "/opt/xt-libsci/11.0.05/cray/74/interlagos/lib/libsci_cray.a")
ELSE()
	# cmake >= 2.8.7 has Hermit.cmake
	SET(CMAKE_SYSTEM_NAME Hermit)
ENDIF()


SET(CMAKE_Fortran_COMPILER ftn)

SET(CMAKE_C_COMPILER cc CACHE FORCE "")
SET(CMAKE_CXX_COMPILER CC CACHE FORCE "")
SET(BUILTIN_LAPACK YES CACHE FORCE "")

# be sure that module xt-libsci is loaded
SET(BUILTIN_BLAS YES CACHE FORCE "")
SET(BUILTIN_MPI YES CACHE FORCE "")

SET(STATIC ON CACHE FORCE "")


# Khabi stuff
#SET(MPI_INCLUDE_PATH $(MPICH_DIR)/include)
#SET(MPI_LIBRARY -L$(MPICH_DIR)/lib)
#set_property(GLOBAL PROPERTY TARGET_SUPPORTS_SHARED_LIBS TRUE)
#set(CMAKE_FIND_LIBRARY_PREFIXES "lib")
#set(CMAKE_FIND_LIBRARY_SUFFIXES ".so" ".a")
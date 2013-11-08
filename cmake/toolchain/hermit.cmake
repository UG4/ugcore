# created by Martin Rupp
# martin.rupp@gcsc.uni-frankfurt.de

# toolchain file for Hermit/XE6 HLRS Stuttgart


# USAGE
# first change programming environment for your compiler of choice
# for Cray compiler:
#    module swap $(module li 2>&1 | awk '/PrgEnv/{print $2}') PrgEnv-cray
# for GCC:
#    module swap $(module li 2>&1 | awk '/PrgEnv/{print $2}') PrgEnv-gnu
# then start cmake
#    cmake -DCMAKE_TOOLCHAIN_FILE=../cmake/toolchain/hermit.cmake ..


# on the cray compiler:
# https://fs.hlrs.de/projects/craydoc/docs_merged/books/S-2179-74/html-S-2179-74/lymwlrwh.html#z862002021malz

# check modules:
# module list
# also: module avail / module load / module swap.

# for use on the CRAY XE6 (Hermit) of HLRS/Karlsruhe
# http://www.hlrs.de/systems/platforms/cray-xe6-hermit/

# see also $CRAY_* environment variables
# $CRAY_UGNI_POST_LINK_OPTS

# to use cmake 2.8.7:
# module load tools/cmake/2.8.7

SET(STATIC_BUILD ON CACHE FORCE "")

SET(HERMIT_CMAKE OFF)
IF(CMAKE_VERSION)
 SET(HERMIT_CMAKE ON)
 IF(${CMAKE_VERSION} STRLESS "2.8.7")
  SET(HERMIT_CMAKE OFF)
  ENDIF()
ENDIF()

IF(HERMIT_CMAKE)
    # cmake >= 2.8.7 has Hermit.cmake
    SET(CMAKE_SYSTEM_NAME Hermit)
ELSE()
    # emulating Hermit.cmake
    # this is needed to get rid of -rdynamic flag...
    SET(CMAKE_SYSTEM_NAME Catamount)
    
    SET(MPI_Fortran_NO_INTERROGATE CMAKE_Fortran_COMPILER)
    SET(MPI_LIBRARY -L$(MPICH_DIR)/lib)
    SET(MPI_EXTRA_LIBRARY -L$(MPICH_DIR)/lib)
    SET(BLA_STATIC ON)
    SET(BLA_VENDOR All)
    SET(BLAS_FIND_QUIETLY ON)
    SET(LAPACK_FIND_QUIETLY ON)
    SET(LAPACK_LIBRARIES "/opt/xt-libsci/11.0.05/cray/74/interlagos/lib/libsci_cray.a")
    
    # Cray is not recognized by Cmake < 2.8.7
    EXECUTE_PROCESS(
	        COMMAND "CC" "-V"
	        OUTPUT_VARIABLE cxx_compiler_string
	        ERROR_VARIABLE cxx_compiler_string
     )
	IF(${cxx_compiler_string} MATCHES ".*\nCray.*")
		SET(CMAKE_CXX_COMPILER_ID "Cray")
		SET(CMAKE_CXX_COMPILER_ID_RUN 1)
	ENDIF()
	
	EXECUTE_PROCESS(
	        COMMAND "cc" "-V"
	        OUTPUT_VARIABLE c_compiler_string
	        ERROR_VARIABLE c_compiler_string
     )
	IF(${c_compiler_string} MATCHES ".*\nCray.*")
		SET(CMAKE_C_COMPILER_ID "Cray")
		SET(CMAKE_C_COMPILER_ID_RUN 1)
	ENDIF()
    
ENDIF()


SET(CMAKE_Fortran_COMPILER ftn)

SET(CMAKE_C_COMPILER cc CACHE FORCE "")
SET(CMAKE_CXX_COMPILER CC CACHE FORCE "")

# be sure that module xt-libsci is loaded
SET(BUILTIN_LAPACK YES CACHE FORCE "")
SET(BUILTIN_BLAS YES CACHE FORCE "")
SET(BUILTIN_MPI YES CACHE FORCE "")



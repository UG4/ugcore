# created by Martin Rupp
# martin.rupp@gcsc.uni-frankfurt.de

# use this file with
# cmake -DCMAKE_TOOLCHAIN_FILE=../toolchain_file_hermit.cmake -DSTATIC=ON ..
# on the cray compiler:
# https://fs.hlrs.de/projects/craydoc/docs_merged/books/S-2179-74/html-S-2179-74/lymwlrwh.html#z862002021malz

# GCC Compiler: For using the gcc, use
# module swap PrgEnv-cray PrgEnv-gnu
# cmake -DHERMIT_USE_GCC=ON -DCMAKE_TOOLCHAIN_FILE=../toolchain_file_hermit.cmake -DSTATIC=ON ..

# check modules:
# module list
# also: module avail / module load / module swap.

# for use on the CRAY XE6 (Hermit) of HLRS/Karlsruhe
# http://www.hlrs.de/systems/platforms/cray-xe6-hermit/

# see also $CRAY_* environment variables
# $CRAY_UGNI_POST_LINK_OPTS


SET(CMAKE_SYSTEM_NAME Catamount)

SET(CMAKE_Fortran_COMPILER ftn)

set(HERMIT_USE_GCC CACHE BOOL "Use GCC on Hermit")
set(CRAY CACHE BOOL "Cray Compiler")

IF(HERMIT_USE_GCC)
 SET(CMAKE_C_COMPILER cc)
 SET(CMAKE_CXX_COMPILER CC)
 SET(BUILTIN_LAPACK)
 SET(BUILTIN_BLAS YES)
 SET(BUILTIN_MPI YES)
 MESSAGE("Info: Using GCC compilers (check output!)")
ELSE(HERMIT_USE_GCC)
 # use the cray compiler
 # check with module list, that PrgEnv-cray is loaded
 SET(CMAKE_C_COMPILER cc)
 SET(CMAKE_CXX_COMPILER CC)
 SET(BUILTIN_LAPACK)
 SET(BUILTIN_BLAS YES)
 SET(BUILTIN_MPI YES)

 SET(CRAY YES)
 MESSAGE("Info: Using the Cray compiler.")
ENDIF(HERMIT_USE_GCC)
SET(STATIC ON)

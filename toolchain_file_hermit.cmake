# created by Martin Rupp
# martin.rupp@gcsc.uni-frankfurt.de

# use this file with
# cmake -DCMAKE_TOOLCHAIN_FILE=../toolchain_file_hermit.cmake -DSTATIC=ON ..

# for use on the CRAY XE6 (Hermit) of HLRS/Karlsruhe
# http://www.hlrs.de/systems/platforms/cray-xe6-hermit/

# https://fs.hlrs.de/projects/craydoc/docs_merged/books/S-2179-74/html-S-2179-74/lymwlrwh.html#z862002021malz

SET(CMAKE_SYSTEM_NAME Catamount)
add_definitions(-DUG_PARALLEL)

# use the cray compiler
# check with module list, that PrgEnv-cray is loaded
SET(CMAKE_C_COMPILER cc)
SET(CMAKE_CXX_COMPILER CC)
SET(CMAKE_Fortran_COMPILER ftn)

# overwrite search for BLAS/LAPACK
add_definitions(-DBLAS_AVAILABLE)
add_definitions(-DLAPACK_AVAILABLE)
SET(LAPACK_FOUND YES)
SET(BLAS_FOUND YES)

# doesn work
SET(STATIC ON)

# created by Martin Rupp
# martin.rupp@gcsc.uni-frankfurt.de

# toolchain file for Hermit/XE6 HLRS Stuttgart with Cray Compiler

# use this file with
# cmake -DCMAKE_TOOLCHAIN_FILE=../toolchain_file_hermit.cmake -DSTATIC=ON ..

# change programming environment
# this does not work 100%, you might have to change it yourself with 
# module swap $(module li 2>&1 | awk '/PrgEnv/{print $2}') PrgEnv-gnu
EXECUTE_PROCESS(COMMAND "module swap $(module li 2>&1 | awk '/PrgEnv/{print $2}') PrgEnv-gnu")

# on the cray compiler:
# https://fs.hlrs.de/projects/craydoc/docs_merged/books/S-2179-74/html-S-2179-74/lymwlrwh.html#z862002021malz

# check modules:
# module list
# also: module avail / module load / module swap.

# for use on the CRAY XE6 (Hermit) of HLRS/Karlsruhe
# http://www.hlrs.de/systems/platforms/cray-xe6-hermit/

# see also $CRAY_* environment variables
# $CRAY_UGNI_POST_LINK_OPTS

SET(CMAKE_SYSTEM_NAME Catamount)
SET(CMAKE_Fortran_COMPILER ftn)

SET(CMAKE_C_COMPILER cc CACHE FORCE "")
SET(CMAKE_CXX_COMPILER CC CACHE FORCE "")
SET(BUILTIN_LAPACK YES CACHE FORCE "")
SET(BUILTIN_BLAS YES CACHE FORCE "")
SET(BUILTIN_MPI YES CACHE FORCE "")

SET(STATIC ON CACHE FORCE "")
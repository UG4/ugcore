# Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
# Author: Martin Rupp
# 
# This file is part of UG4.
# 
# UG4 is free software: you can redistribute it and/or modify it under the
# terms of the GNU Lesser General Public License version 3 (as published by the
# Free Software Foundation) with the following additional attribution
# requirements (according to LGPL/GPL v3 §7):
# 
# (1) The following notice must be displayed in the Appropriate Legal Notices
# of covered and combined works: "Based on UG4 (www.ug4.org/license)".
# 
# (2) The following notice must be displayed at a prominent place in the
# terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
# 
# (3) The following bibliography is recommended for citation and must be
# preserved in all covered files:
# "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
#   parallel geometric multigrid solver on hierarchically distributed grids.
#   Computing and visualization in science 16, 4 (2013), 151-164"
# "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
#   flexible software system for simulating pde based models on high performance
#   computers. Computing and visualization in science 16, 4 (2013), 165-179"
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Lesser General Public License for more details.

# included from build_algebra.cmake
########################################

# if the search for LAPACK/BLAS is not working, you can use
# a) USER_LAPACK_LIBRARIES, (optionally  USER_LAPACK_INCLUDE_DIR ) and USER_BLAS_LIBRARIES (optionally USER_BLAS_INCLUDE_DIR )
# b) BUILTIN_LAPACK and BUILTIN_BLAS (this is used e.g. on the Hermit cluster)
#
# Some common paths for LAPACK/BLAS are:
# /usr/local/lib /usr/lib /usr/local/lib64 /usr/lib64 /bgsys/local/lib
# Please mail martin.rupp@gcsc.uni-frankfurt.de if you have problems setting up
# OR if you solved your problem.

if(STATIC_BUILD)
	set(BLA_STATIC ON)
endif()

# If LAPACK is enabled, try to find the library
if(LAPACK)
    # try find using cmake-native find_package (requires fortran)
	find_package(LAPACK QUIET)

    # if not found, try to used ug4-shipped c-version (no fortran requirement)
	if(NOT LAPACK_FOUND)
    	find_package(MYLAPACK QUIET)
	endif(NOT LAPACK_FOUND)

    # now display search status:
    # a) User defined libraries are to used. (Show found libraries for info)
	if(USER_LAPACK_LIBRARIES)
	    message(STATUS "Info: Using USER-defined LAPACK (Include: ${USER_LAPACK_INCLUDE_DIR}, Lib: ${USER_LAPACK_LIBRARIES})")
	    if(LAPACK_FOUND)
        message(STATUS "Note: CMake also found (unused) LAPACK (Includes: ${LAPACK_INCLUDE_DIR}, Lib: ${LAPACK_LIBRARIES})")
	    endif(LAPACK_FOUND)
		include_directories (${USER_LAPACK_INCLUDE_DIR})
		set(linkLibraries ${linkLibraries} ${USER_LAPACK_LIBRARIES})
		add_definitions(-DLAPACK_AVAILABLE)	
	# b) Using build-in LAPACK (i.e. nothing needed, added by e.g. compiler)
	elseif(BUILTIN_LAPACK)
		message(STATUS "Info: Using Builtin LAPACK")
		add_definitions(-DLAPACK_AVAILABLE)
	# c) Lapack has been found
	elseif(LAPACK_FOUND)
		message(STATUS "Info: Using LAPACK (Include: ${LAPACK_INCLUDE_DIR}, Lib: ${LAPACK_LIBRARIES})")
		include_directories (${LAPACK_INCLUDE_DIR})
		set(linkLibraries ${linkLibraries} ${LAPACK_LIBRARIES})
		add_definitions(-DLAPACK_AVAILABLE)	
	# d) Not found
	else()	
		message(STATUS "WARNING: No LAPACK package found. Not using LAPACK.")
		message(STATUS "         If you have a builtin LAPACK, use -DBUILTIN_LAPACK")
		message(STATUS "         To supply your own LAPACK libraries, use -DUSER_LAPACK_LIBRARIES and (optionally) -DUSER_LAPACK_INCLUDES:")
		message(STATUS "         e.g.: cmake -DUSER_LAPACK_LIBRARIES=/usr/lib64/liblapack.so -DUSER_LAPACK_INCLUDES=/includepath/ ..")
	endif()
else(LAPACK)
    # e) Lapack not wanted
	message(STATUS "Info: Not using LAPACK, use -DLAPACK=ON to enable.")
endif(LAPACK)

# If BLAS is enabled, try to find the library
if(BLAS)
    # try find using cmake-native find_package (requires fortran)
	find_package(BLAS QUIET)

    # if not found, try to used ug4-shipped c-version (no fortran requirement)
	if(NOT BLAS_FOUND)
    	find_package(MYBLAS QUIET)
	endif(NOT BLAS_FOUND)
	
	# this is to fix something. todo: add comment
	if(BLAS_goto2_LIBRARY)
		message(STATUS "Info: GotoBLAS2 found: (${BLAS_goto2_LIBRARY}). Adding -lgfortran.")
		SET(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -lgfortran")
	endif()	
	
    # now display search status:
    # a) User defined libraries are to used. (Show found libraries for info)
	if(USER_BLAS_LIBRARIES)
		message(STATUS "Info: Using USER-defined BLAS (Include: ${USER_BLAS_INCLUDE_DIR}, Lib: ${USER_BLAS_LIBRARIES})")
    	if(BLAS_FOUND)
	    message(STATUS "Note: CMake also found (unused) BLAS (Includes: ${BLAS_INCLUDE_DIR}, Lib: ${BLAS_LIBRARIES})")
	    endif(BLAS_FOUND)
		include_directories (${USER_BLAS_INCLUDE_DIR})
		set(linkLibraries ${linkLibraries} ${USER_BLAS_LIBRARIES})
		add_definitions(-DBLAS_AVAILABLE)
	# b) Using build-in BLAS (i.e. nothing needed, added by e.g. compiler)
	elseif(BUILTIN_BLAS)
		message(STATUS "Info: Using Builtin BLAS")
		add_definitions(-DBLAS_AVAILABLE)
	# c) Blas has been found
	elseif(BLAS_FOUND)
		message(STATUS "Info: Using BLAS (Include: ${BLAS_INCLUDE_DIR}, Lib: ${BLAS_LIBRARIES})")
		include_directories (${BLAS_INCLUDE_DIR})
		set(linkLibraries ${linkLibraries} ${BLAS_LIBRARIES})
		add_definitions(-DBLAS_AVAILABLE)
	# d) Not found
	else()	
		message(STATUS "WARNING: No BLAS package found. Not using BLAS.")
		message(STATUS "         If you have a builtin BLAS, use -DBUILTIN_BLAS")
		message(STATUS "         To supply your own BLAS libraries, use -DUSER_BLAS_LIBRARIES and (optionally) -DUSER_BLAS_INCLUDES:")
		message(STATUS "         e.g.: cmake -DUSER_BLAS_LIBRARIES=/usr/lib64/liblablas.so -DUSER_BLAS_INCLUDES=/includepath/ ..")
	endif()
else(BLAS)
    # e) Blas not wanted
	message(STATUS "Info: Not using BLAS, use -DBLAS=ON to enable.")
endif(BLAS)

################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_algebra
################################################

cmake_minimum_required (VERSION 2.6)

project (P_LIB_ALGEBRA)


if(WITH_ALGEBRA)

#######################
# Arne's Algebra

find_package(Boost)

if(Boost_FOUND)
	message("-- Boost - Library found at '${Boost_INCLUDE_DIRS}'. \n   Compiling WITH Boost.")

	# add include path	
	include_directories (../solver/ ${Boost_INCLUDE_DIRS})

else(Boost_FOUND)
	message("-- WARNING: Boost - Libraray NOT found. Compiling WITHOUT Boost.")
endif(Boost_FOUND)


########################
# search hyprelib

find_path (HYPRELIB_LIB_DIR libHYPRE.a 
		PATHS ENV PATH
		PATH_SUFFIXES src/hypre/lib/ )

if(HYPRELIB_LIB_DIR)
	message("-- HYPRE - Library found at '${HYPRELIB_LIB_DIR}'.\n   Compiling WITH Hypre.")

	# add include path
	include_directories (${HYPRELIB_LIB_DIR}/../include/)

	# add lib path
	link_directories(${HYPRELIB_LIB_DIR})

	# let all executables link against libHYPRE.a 
	link_libraries (${HYPRELIB_LIB_DIR}/libHYPRE.a)

	# add include directories
	include_directories(${P_UGBASE_SOURCE_DIR})

	set(srcHypre_Algebra	 hypre_algebra/matrix.cpp
		     		 hypre_algebra/vector.cpp
	       			 hypre_algebra/linearsolver.cpp)

	add_library(algebra ${srcHypre_Algebra} )
	add_dependencies(algebra common)

else(HYPRELIB_LIB_DIR)
	message("-- WARNING: HYPRE - Libraray NOT found. Compiling WITHOUT Hypre.")
endif(HYPRELIB_LIB_DIR)


##########################
# setup library

endif(WITH_ALGEBRA)


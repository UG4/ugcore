################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_algebra
################################################

cmake_minimum_required (VERSION 2.6)

project (P_LIB_ALGEBRA)


if(WITH_ALGEBRA)

if(Boost_FOUND)
	set(srcHypre_Algebra	 arne_algebra/arnevector.cpp
				 arne_algebra/arnematrix.cpp
				 arne_algebra/arnelinearsolver.cpp)
endif(Boost_FOUND)

if(HYPRELIB_LIB_DIR)

	set(srcHypre_Algebra	 ${srcHypre_Algebra}
				 hypre_algebra/hyprematrix.cpp
		     		 hypre_algebra/hyprevector.cpp
	       			 hypre_algebra/hyprelinearsolver.cpp)
endif(HYPRELIB_LIB_DIR)

##########################
# setup library

add_library(algebra ${srcHypre_Algebra} )
add_dependencies(algebra common lib_algebra.h)

endif(WITH_ALGEBRA)


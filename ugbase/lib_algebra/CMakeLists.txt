################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_algebra
################################################

cmake_minimum_required (VERSION 2.6)

project (P_LIB_ALGEBRA)

include_directories(${P_UGBASE_SOURCE_DIR})

if(NOT NO_ALGEBRA)
	set(src_Algebra)

	##########################
	# ublas algebra
	##########################
	if(Boost_FOUND)
		set(src_Algebra	 ublas_algebra/ublas_vector.cpp
				 ublas_algebra/ublas_matrix.cpp
				 ublas_algebra/ublas_linearsolver.cpp)
	endif(Boost_FOUND)

	##########################
	# hypre algebra
	##########################
	if(HYPRELIB_DIR)
		set(src_Algebra	 ${src_Algebra}
				 hypre_algebra/hyprematrix.cpp
		     		 hypre_algebra/hyprevector.cpp
	       		 	 hypre_algebra/hyprelinearsolver.cpp)
	endif(HYPRELIB_DIR)

	##########################
	# martin algebra
	##########################
#	set(src_Algebra	 ${src_Algebra}
#			 martin_algebra/misc.cpp)


	##########################

	set(src_Algebra	 ${src_Algebra}
		operator/preconditioner/amg/amg_coarsening.cpp )

	##########################
	# setup library
        EXPORTSOURCES(lib_algebra ${src_Algebra})

if(NOT UG_SHARED)
	add_library(algebra ${src_Algebra} )
	add_dependencies(algebra common lib_algebra.h)
endif(NOT UG_SHARED)

endif(NOT NO_ALGEBRA)


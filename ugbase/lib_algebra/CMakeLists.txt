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
	# cpu algebra
	##########################


	##########################

	set(src_Algebra	 ${src_Algebra}
		operator/preconditioner/amg/amg_coarsening.cpp )

if(NOT NO_PCL)
	# add parallelization
	set(src_Algebra	 ${src_Algebra}
						operator/linear_solver/feti.cpp
						parallelization/parallelization_util.cpp)
endif(NOT NO_PCL)

	##########################
	# setup library
        EXPORTSOURCES(lib_algebra ${src_Algebra})

if(NOT UG_SHARED)
	add_library(algebra ${src_Algebra} )
	#add_dependencies(algebra common lib_algebra.h)
endif(NOT UG_SHARED)

endif(NOT NO_ALGEBRA)


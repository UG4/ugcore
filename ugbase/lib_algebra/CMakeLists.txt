################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_algebra
################################################

cmake_minimum_required (VERSION 2.6)

project (P_LIB_ALGEBRA)

include_directories(${P_UGBASE_SOURCE_DIR})

set(src_Algebra)

##########################
# cpu algebra
##########################


##########################

set(src_Algebra	 ${src_Algebra}
	algebra_type.cpp
	common/connection_viewer_output.cpp
	common/connection_viewer_input.cpp
	 )
add_subdirectory(common/matrixio)

if(PARALLEL)
	# add parallelization
	set(src_Algebra	 ${src_Algebra}
						operator/linear_solver/feti.cpp
						parallelization/parallelization_util.cpp
						parallelization/parallel_coloring.cpp						
						parallelization/global_layout.cpp
						parallelization/parallel_index_layout.cpp
						parallelization/parallel_nodes.cpp	
						 )
endif(PARALLEL)

##########################
# setup library
if(BUILD_ONE_LIB)
    EXPORTSOURCES(ugbase/lib_algebra ${src_Algebra})
else(BUILD_ONE_LIB)
	add_library(algebra ${src_Algebra} )
endif(BUILD_ONE_LIB)



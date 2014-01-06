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
	small_algebra/solve_deficit.cpp
	operator/preconditioner/line_smoothers.cpp
	algebra_common/permutation_util.cpp
	)
add_subdirectory(common/matrixio)

if(PARALLEL)
	# add parallelization
	set(src_Algebra	 ${src_Algebra}
						operator/preconditioner/schur.cpp
	    				operator/linear_solver/feti.cpp
						parallelization/parallelization_util.cpp
						parallelization/parallel_coloring.cpp						
						parallelization/global_layout.cpp
						parallelization/parallel_index_layout.cpp
						parallelization/parallel_nodes.cpp	
						parallelization/algebra_layouts.cpp						
						 )
endif(PARALLEL)

if(CUDA_FOUND)
	set(src_Algebra ${src_Algebra} 
			 	gpu_algebra/cuda/cuda_error.cpp
				gpu_algebra/cuda/cuda_helper.cpp                
                gpu_algebra/cuda/helper_cuda.cpp
                gpu_algebra/cuda/common_cuda.cu
				)
endif(CUDA_FOUND)

##########################
# setup library
if(BUILD_ONE_LIB)
  EXPORTSOURCES(ugbase/lib_algebra ${src_Algebra})   
else(BUILD_ONE_LIB)
    CUDA_ADD_LIBRARY( lib_algebra SHARED ${src_Algebra}  )
   #add_library(algebra_CUDA ${src_Algebra})
endif(BUILD_ONE_LIB)

# Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
# Author: Andreas Vogel
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

cmake_minimum_required (VERSION 2.6)

project (P_LIB_ALGEBRA)

include_directories(${P_UGBASE_SOURCE_DIR})

set(src_Algebra)

##########################
# cpu algebra
##########################


##########################

set(src_Algebra	 ${src_Algebra}
    debug_ids.cpp
	algebra_type.cpp
	common/connection_viewer_output.cpp
	common/connection_viewer_input.cpp
	small_algebra/solve_deficit.cpp
	operator/preconditioner/line_smoothers.cpp
	operator/linear_solver/analyzing_solver.cpp
	algebra_common/permutation_util.cpp
	operator/preconditioner/schur/schur.cpp
	)
	
add_subdirectory(common/matrixio)

if(PARALLEL)
	# add parallelization
	set(src_Algebra	 ${src_Algebra}
						operator/preconditioner/schur/schur_precond.cpp
						operator/preconditioner/schur/schur_complement_operator.cpp
						operator/preconditioner/schur/parallel_progress.cpp
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
				gpu_algebra/cuda/cuda_manager.cpp                
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

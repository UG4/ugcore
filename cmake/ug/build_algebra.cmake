# included from ug_includes.cmake
########################################
# buildAlgebra
if(buildAlgebra)
	add_definitions(-DUG_ALGEBRA)
	
	include(${UG_ROOT_PATH}/cmake/ug/lapack_blas.cmake)
	
endif(buildAlgebra)
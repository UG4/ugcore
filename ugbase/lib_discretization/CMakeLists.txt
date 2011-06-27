################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include("../../ug_cmake_includes.txt")

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcDiscretization	common/function_group.cpp
						common/subset_group.cpp
						common/groups_util.cpp

						dof_manager/function_pattern.cpp
						dof_manager/cuthill_mckee.cpp
						dof_manager/p1conform/p1conform.cpp

						reference_element/reference_element.cpp
						
						quadrature/gauss_quad/gauss_quad_edge.cpp
						quadrature/gauss_quad/gauss_quad_triangle.cpp
						quadrature/gauss_quad/gauss_quad_quadrilateral.cpp
						quadrature/gauss_quad/gauss_quad_tetrahedron.cpp
						quadrature/gauss_quad/gauss_quad_prism.cpp
						quadrature/gauss_quad/gauss_quad_pyramid.cpp
						quadrature/gauss_quad/gauss_quad_hexahedron.cpp

						local_shape_function_set/local_shape_function_set_provider.cpp
						local_shape_function_set/lagrange/lagrangep1.cpp)

# add parallelization
if(NOT NO_MPI)
	set(srcParallelization	parallelization/parallelization_util.cpp)
else(NOT NO_MPI)
	set(srcParallelization )
endif(NOT NO_MPI)

set(Sources         ${srcParallelization}
					${srcDiscretization})

EXPORTSOURCES(lib_discretization ${Sources})

if(NOT UG_SHARED)
add_library(discretization STATIC ${Sources})
#add_dependencies(discretization algebra grid common)
endif(NOT UG_SHARED)

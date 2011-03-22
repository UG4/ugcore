################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcDiscretization	common/function_group.cpp
						common/subset_group.cpp
						common/groups_util.cpp

						dof_manager/function_pattern.cpp
						dof_manager/cuthill_mckee.cpp
						dof_manager/p1conform/p1conform.cpp

						reference_element/reference_element.cpp
						
						quadrature/gauss_quadratures_ReferenceEdge.cpp
						quadrature/gauss_quadratures_ReferenceTriangle.cpp
						quadrature/gauss_quadratures_ReferenceQuadrilateral.cpp
						quadrature/gauss_quadratures_ReferenceTetrahedron.cpp
						quadrature/gauss_quadratures_ReferencePrism.cpp
						quadrature/gauss_quadratures_ReferencePyramid.cpp
						quadrature/gauss_quadratures_ReferenceHexahedron.cpp

						local_shape_function_set/local_shape_function_set_provider.cpp
						local_shape_function_set/local_dof_pattern.cpp
						local_shape_function_set/LagrangeP1/lagrangep1.cpp)


EXPORTSOURCES(lib_discretization ${srcDiscretization})

if(NOT UG_SHARED)
add_library(discretization STATIC ${srcDiscretization})
#add_dependencies(discretization algebra grid common)
endif(NOT UG_SHARED)

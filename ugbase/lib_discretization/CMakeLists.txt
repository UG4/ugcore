################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcDiscretization		dof_manager/general_dof_manager/general_dof_manager.cpp
				dof_manager/p1conform_dof_manager/p1conform_dof_manager.cpp
				quadrature/gauss_quadratures_ReferenceTriangle.cpp
				quadrature/gauss_quadratures_ReferenceQuadrilateral.cpp
				local_shape_function_set/local_shape_function_set_factory.cpp
				local_shape_function_set/local_dof_pattern.cpp
				local_shape_function_set/LagrangeP1/p1conform.cpp)

				#quadrature/gauss_quadratures_ReferenceVertex.cpp
				#quadrature/gauss_quadratures_ReferenceTetrahedron.cpp
				#quadrature/gauss_quadratures_ReferencePrism.cpp
				#quadrature/gauss_quadratures_ReferencePyramid.cpp
				#quadrature/gauss_quadratures_ReferenceHexahedron.cpp

				#vtkoutput.cpp
			     	#newton.cpp)
				#differentialoperator.cpp
				#elementdata.cpp
				#equation.cpp

add_library(discretization 	${srcDiscretization})

add_dependencies(discretization algebra grid common)

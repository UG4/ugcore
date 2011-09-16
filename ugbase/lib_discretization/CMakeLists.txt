################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include("../../ug_cmake_includes.txt")

set(srcDiscretization	common/function_group.cpp
						common/subset_group.cpp
						common/groups_util.cpp

						dof_manager/function_pattern.cpp
						dof_manager/cuthill_mckee.cpp
						dof_manager/lexorder.cpp
						dof_manager/p1conform/p1conform.cpp
						dof_manager/conform/conform.cpp

						reference_element/reference_element.cpp
			            reference_element/reference_element_init.cpp
			            reference_element/reference_element_mapping.cpp
						
						quadrature/gauss_quad/gauss_quad_edge.cpp
						quadrature/gauss_quad/gauss_quad_triangle.cpp
						quadrature/gauss_quad/gauss_quad_quadrilateral.cpp
						quadrature/gauss_quad/gauss_quad_tetrahedron.cpp
						quadrature/gauss_quad/gauss_quad_prism.cpp
						quadrature/gauss_quad/gauss_quad_pyramid.cpp
						quadrature/gauss_quad/gauss_quad_hexahedron.cpp

						local_finite_element/local_finite_element_id.cpp
						local_finite_element/local_shape_function_set.cpp
						local_finite_element/local_dof_set.cpp
						local_finite_element/lagrange/lagrangep1.cpp
						local_finite_element/lagrange/lagrange.cpp
						
						operator/linear_operator/multi_grid_solver/mg_solver_util.cpp
						
						spatial_discretization/elem_disc/elem_disc_interface.cpp
						spatial_discretization/disc_util/finite_volume_geometry.cpp
						spatial_discretization/disc_util/hanging_finite_volume_geometry.cpp
						spatial_discretization/ip_data/data_evaluator.cpp
						
						spatial_discretization/elem_disc/constant_equation/constant_equation_fv1.cpp
						spatial_discretization/elem_disc/convection_diffusion/convection_diffusion.cpp
						spatial_discretization/elem_disc/density_driven_flow/density_driven_flow_fv1.cpp
						spatial_discretization/elem_disc/inner_boundary/inner_boundary_fv1.cpp
						spatial_discretization/elem_disc/linear_elasticity/linear_elasticity_fe1.cpp
						spatial_discretization/elem_disc/neumann_boundary/neumann_boundary_fv1.cpp
						spatial_discretization/elem_disc/nonlinear_elasticity/nonlinear_elasticity_fe1.cpp
						spatial_discretization/elem_disc/thermohaline_flow/thermohaline_flow_fv1.cpp
						spatial_discretization/elem_disc/time_neumann_boundary/time_neumann_boundary.cpp)
						

# add parallelization
if(PARALLEL)
	set(srcParallelization	parallelization/parallelization_util.cpp)
else(PARALLEL)
	set(srcParallelization )
endif(PARALLEL)

set(Sources         ${srcParallelization}
					${srcDiscretization})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(lib_discretization ${Sources})
else(BUILD_ONE_LIB)
	add_library(discretization STATIC ${Sources})
endif(BUILD_ONE_LIB)

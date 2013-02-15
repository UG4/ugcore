################################################
# created by Andreas Vogel
# andreas.vogel@gcsc.uni-frankfurt.de
#
# lib_discretization
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include("../../ug_cmake_includes.txt")

set(srcDiscretization	domain_util.cpp

                        common/function_group.cpp
						common/subset_group.cpp
						common/groups_util.cpp
						common/subset_util.cpp

						dof_manager/function_pattern.cpp
						dof_manager/cuthill_mckee.cpp
						dof_manager/lexorder.cpp
						dof_manager/dof_distribution_info.cpp
						dof_manager/mg_dof_distribution.cpp
						dof_manager/level_dof_distribution.cpp
						dof_manager/surface_dof_distribution.cpp
						dof_manager/managing_dof_distribution.cpp

                        function_spaces/approximation_space.cpp
                        function_spaces/dof_position_util.cpp
                        function_spaces/grid_function.cpp

                        io/vtkoutput.cpp

						reference_element/reference_element.cpp
			            reference_element/reference_mapping_provider.cpp
						
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
					
						spatial_disc/subset_assemble_util.cpp
						spatial_disc/elem_disc/elem_disc_interface.cpp
						spatial_disc/disc_util/fvho_geom.cpp
						spatial_disc/disc_util/fv1_geom.cpp
						spatial_disc/disc_util/fvcr_geom.cpp
						spatial_disc/disc_util/hfv1_geom.cpp
						spatial_disc/user_data/data_evaluator.cpp
						
						spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.cpp
						
						spatial_disc/elem_disc/inner_boundary/inner_boundary_fv1.cpp
						spatial_disc/elem_disc/neumann_boundary/neumann_boundary.cpp)
						
						
						

# add parallelization
if(PARALLEL)
	set(srcParallelization	parallelization/parallelization_util.cpp)
else(PARALLEL)
	set(srcParallelization )
endif(PARALLEL)

set(Sources         ${srcParallelization}
					${srcDiscretization})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/lib_disc ${Sources})
else(BUILD_ONE_LIB)
	add_library(disc STATIC ${Sources})
endif(BUILD_ONE_LIB)

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

cmake_minimum_required(VERSION 2.6)

project(P_LIB_DISCRETIZATION)

include("../../cmake/ug_includes.cmake")

set(srcDiscretization	domain.cpp
						domain_util.cpp

                        common/function_group.cpp
						common/groups_util.cpp

						dof_manager/function_pattern.cpp
						dof_manager/orientation.cpp
						dof_manager/dof_count.cpp
						dof_manager/dof_index_storage.cpp
						dof_manager/dof_distribution_info.cpp
						dof_manager/dof_distribution.cpp
						dof_manager/ordering/cuthill_mckee.cpp
						dof_manager/ordering/lexorder.cpp
						dof_manager/ordering/downwindorder.cpp

                        function_spaces/approximation_space.cpp
                        function_spaces/dof_position_util.cpp
                        function_spaces/grid_function.cpp
                        function_spaces/adaption_surface_grid_function.cpp
                        function_spaces/local_transfer_interface.cpp

                        io/vtkoutput.cpp

						reference_element/reference_element.cpp
			            reference_element/reference_mapping_provider.cpp
						
						quadrature/quad_test.cpp
						quadrature/quadrature_provider.cpp
						quadrature/gauss_legendre/gauss_legendre.cpp
						quadrature/newton_cotes/newton_cotes.cpp
						quadrature/gauss_jacobi/gauss_jacobi20.cpp
						quadrature/gauss_jacobi/gauss_jacobi10.cpp
						quadrature/gauss_tensor_prod/gauss_tensor_prod.cpp
						quadrature/gauss/gauss_quad_vertex.cpp
						quadrature/gauss/gauss_quad_edge.cpp
						quadrature/gauss/gauss_quad_triangle.cpp
						quadrature/gauss/gauss_quad_quadrilateral.cpp
						quadrature/gauss/gauss_quad_tetrahedron.cpp
						quadrature/gauss/gauss_quad_prism.cpp
						quadrature/gauss/gauss_quad_pyramid.cpp
						quadrature/gauss/gauss_quad_hexahedron.cpp
						quadrature/gauss/gauss_quad_octahedron.cpp

						local_finite_element/lagrange/lagrange_local_dof.cpp
						local_finite_element/lagrange/lagrangep1.cpp
						local_finite_element/lagrange/lagrange.cpp
						local_finite_element/local_finite_element_id.cpp
						local_finite_element/local_finite_element_provider.cpp
						local_finite_element/local_dof_set.cpp
						local_finite_element/mini/mini.cpp
						
						operator/linear_operator/multi_grid_solver/mg_solver.cpp
						
						spatial_disc/subset_assemble_util.cpp
						spatial_disc/elem_disc/elem_disc_interface.cpp
						spatial_disc/disc_util/fe_geom.cpp
						spatial_disc/disc_util/fvho_geom.cpp
						spatial_disc/disc_util/fv1_geom.cpp
                        spatial_disc/disc_util/fv1FT_geom.cpp
            spatial_disc/disc_util/fv1Cut_geom.cpp
						spatial_disc/disc_util/fvcr_geom.cpp
						spatial_disc/disc_util/hfv1_geom.cpp
						spatial_disc/disc_util/hfvcr_geom.cpp
						spatial_disc/user_data/data_evaluator.cpp
						spatial_disc/user_data/data_export.cpp
						
						spatial_disc/constraints/continuity_constraints/p1_continuity_constraints.cpp
						
						spatial_disc/elem_disc/inner_boundary/inner_boundary_fv1.cpp
						spatial_disc/elem_disc/neumann_boundary/neumann_boundary_base.cpp
						spatial_disc/elem_disc/neumann_boundary/fv1/neumann_boundary_fv1.cpp
						spatial_disc/elem_disc/neumann_boundary/fe/neumann_boundary_fe.cpp
						spatial_disc/elem_disc/neumann_boundary/fv/neumann_boundary_fv.cpp)
						
						
						

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

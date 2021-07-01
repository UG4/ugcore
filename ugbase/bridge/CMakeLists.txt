# Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
# Author: Sebastian Reiter
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

####
# ug_bridge Library
####

project(P_UGBRIDGE)

include("../../cmake/ug_includes.cmake")

set(srcUGBridge    	bridge.cpp
					misc_bridges/test_bridge.cpp
					misc_bridges/profiler_bridge.cpp
					misc_bridges/util_bridge.cpp
					misc_bridges/misc_bridge.cpp
					misc_bridges/raster_bridge.cpp
					misc_bridges/vec_math_bridge.cpp
					misc_bridges/orthopoly_bridge.cpp)

set(srcUGBridge     ${srcUGBridge}
					misc_bridges/pcl_bridge.cpp)

if(buildGrid)
	set(srcUGBridge     ${srcUGBridge}
						grid_bridges/balancing_bridge.cpp
						grid_bridges/debug_bridge.cpp
						grid_bridges/file_io_bridge.cpp
						grid_bridges/grid_bridge.cpp
						grid_bridges/grid_bridges.cpp
						grid_bridges/layers_bridge.cpp
						grid_bridges/misc_bridge.cpp
						grid_bridges/refinement_bridge.cpp
						grid_bridges/selector_bridge.cpp
						grid_bridges/subset_handler_bridge.cpp)
endif(buildGrid)

# if algebra is enabled, we have to add additional src-files and link with the
# discretization module.
if(buildAlgebra)
	# add source files that depend on the algebra and discretization system here.
	set(srcUGBridge ${srcUGBridge}
					domain_bridges/domain_bridge.cpp
					domain_bridges/load_balancing_bridge.cpp
					domain_bridges/periodic_boundary_bridge.cpp
					domain_bridges/ray_tracing_bridge.cpp
					domain_bridges/refinement_bridge.cpp
					domain_bridges/selection_bridge.cpp
					domain_bridges/transform_bridge.cpp

					algebra_bridges/common_bridge.cpp
					algebra_bridges/preconditioner_bridge.cpp
					algebra_bridges/schur_bridge.cpp
					algebra_bridges/ordering_bridge.cpp
					algebra_bridges/solver_bridge.cpp
					algebra_bridges/eigensolver_bridge.cpp
					algebra_bridges/domain_dependent_preconditioner_bridge.cpp
					#algebra_bridges/constrained_linear_iterator_bridge.cpp
					algebra_bridges/pilut_bridge.cpp
					algebra_bridges/obstacle_bridge.cpp
					algebra_bridges/restart_bridge.cpp
					
					disc_bridges/domain_disc_bridge.cpp
					disc_bridges/interpolate_bridge.cpp
					disc_bridges/evaluate_at_position_bridge.cpp					
					disc_bridges/max_error_bridge.cpp
					disc_bridges/grid_function_bridge.cpp
					disc_bridges/ordering_bridge.cpp
					disc_bridges/integrate_bridge.cpp
					disc_bridges/adaptive_tools_bridge.cpp
					disc_bridges/common_bridge.cpp
					disc_bridges/algebra_bridge.cpp
					disc_bridges/elem_discs_bridge.cpp
					disc_bridges/constraints_bridge.cpp
					disc_bridges/multigrid_bridge.cpp
					disc_bridges/output_bridge.cpp
					disc_bridges/finite_volume_bridge.cpp
					disc_bridges/user_data_bridge.cpp
					disc_bridges/manifold_util_bridge.cpp
					disc_bridges/reference_mapping_test_bridge.cpp)
else(buildAlgebra)

endif(buildAlgebra)

set(Sources ${srcUGBridge})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/bridge  ${Sources})
else(BUILD_ONE_LIB)
	add_library(bridge  ${Sources})
endif(BUILD_ONE_LIB)

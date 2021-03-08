# Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

################################################
# lib_grid
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_GRID)

include("../../cmake/ug_includes.cmake")

set(srcGrid		grid/grid.cpp
				grid/grid_base_objects.cpp
				grid/grid_connection_managment.cpp
				grid/grid_object_collection.cpp
				grid/grid_util.cpp
				grid/neighborhood.cpp
				grid/neighborhood_util.cpp)
				
set(srcAlgorithms	algorithms/debug_util.cpp
					algorithms/element_side_util.cpp
					algorithms/field_util.cpp
					algorithms/grid_statistics.cpp
					algorithms/heightfield_util.cpp
					algorithms/hexahedron_util.cpp
					algorithms/subset_util.cpp
					algorithms/subset_dim_util.cpp
					algorithms/selection_util.cpp
					algorithms/serialization.cpp
					algorithms/orientation_util.cpp
					algorithms/polychain_util.cpp
					algorithms/problem_detection_util.cpp
					algorithms/quadrilateral_util.cpp
					algorithms/quality_util.cpp
					algorithms/raster_layer_util.cpp
					algorithms/ray_element_intersection_util.cpp
					algorithms/subset_color_util.cpp
					algorithms/remeshing/delaunay_info.cpp
					algorithms/remeshing/delaunay_triangulation.cpp
					algorithms/remeshing/edge_length_adjustment.cpp
					algorithms/remeshing/edge_length_adjustment_extended.cpp
					algorithms/remeshing/simple_grid.cpp
					algorithms/remeshing/grid_adaption.cpp
					algorithms/geom_obj_util/edge_util.cpp
					algorithms/geom_obj_util/face_util.cpp
					algorithms/geom_obj_util/misc_util.cpp
					algorithms/geom_obj_util/vertex_util.cpp
					algorithms/geom_obj_util/volume_util.cpp
					algorithms/grid_generation/horizontal_layers_mesher.cpp
					algorithms/grid_generation/icosahedron.cpp
					algorithms/grid_generation/tetrahedralization.cpp
					algorithms/grid_generation/triangle_fill.cpp
					algorithms/grid_generation/triangle_fill_sweep_line.cpp
					algorithms/extrusion/extrude.cpp
					algorithms/extrusion/cylinder_extrusion.cpp
					algorithms/extrusion/expand_layers.cpp
					algorithms/projections/overlying_subset_finder.hpp
					algorithms/projections/z_ray_tracer.hpp
					refinement/projectors/projection_handler.cpp
					refinement/projectors/smooth_projector.cpp
					refinement/projectors/subdivision_projector.cpp
					refinement/projectors/neurite_projector.cpp
					refinement/projectors/elliptic_cylinder_projector.cpp
					refinement/ref_mark_adjusters/local_mark_adjuster.cpp
					refinement/ref_mark_adjusters/horizontal_anisotropy_adjuster.cpp
					refinement/ref_mark_adjusters/shadow_copy_adjuster.cpp
					refinement/ref_mark_adjusters/mg_hnode_adjuster.cpp
					refinement/ref_mark_adjusters/std_hnode_adjuster.cpp
					refinement/adaptive_regular_mg_refiner.cpp
					refinement/fractured_media_refiner.cpp
					refinement/hanging_node_refiner_base.cpp
					refinement/hanging_node_refiner_grid.cpp
					refinement/hanging_node_refiner_multi_grid.cpp
					refinement/refiner_interface.cpp
					refinement/regular_refinement.cpp
					refinement/global_multi_grid_refiner.cpp
					refinement/global_subdivision_multi_grid_refiner.cpp
					refinement/global_fractured_media_refiner.cpp
					algorithms/subdivision/subdivision_loop.cpp
					algorithms/subdivision/subdivision_rules_piecewise_loop.cpp
					algorithms/subdivision/subdivision_volumes.cpp
					algorithms/tkd/tkd_info.cpp
					algorithms/tkd/tkd_util.cpp
					algorithms/unit_tests/check_associated_elements.cpp)
					
set(srcFileIO	file_io/file_io_2df.cpp
    			file_io/file_io_art.cpp
				file_io/file_io_asc.cpp
				file_io/file_io_dump.cpp
				file_io/file_io_lgb.cpp
				file_io/file_io_lgm.cpp
				file_io/file_io_ng.cpp
				file_io/file_io_obj.cpp
				file_io/file_io_tetgen.cpp
				file_io/file_io_tikz.cpp
				file_io/file_io_txt.cpp
				file_io/file_io_ug.cpp
				file_io/file_io_ugx.cpp
				file_io/file_io_ncdf.cpp
				file_io/file_io_msh.cpp
				file_io/file_io_stl.cpp
				file_io/file_io_vtu.cpp
				file_io/file_io_swc.cpp
				file_io/file_io.cpp)
								
set(srcLibGrid	common_attachments.cpp
                tools/bool_marker.cpp
				tools/subset_handler_interface.cpp
				tools/subset_handler_grid.cpp
				tools/subset_handler_multi_grid.cpp
				tools/selector_interface.cpp
				tools/selector_grid.cpp
				tools/selector_multi_grid.cpp
				tools/marker_points.cpp
				tools/surface_view.cpp
				tools/partition_map.cpp
				tools/periodic_boundary_manager.cpp
				tools/grid_level.cpp
				tools/subset_group.cpp
				grid_objects/grid_objects_1d.cpp
				grid_objects/grid_objects_2d.cpp
				grid_objects/grid_objects_3d.cpp
				grid_objects/tetrahedron_rules.cpp
				grid_objects/octahedron_rules.cpp
				grid_objects/pyramid_rules.cpp
				grid_objects/prism_rules.cpp
				grid_objects/hexahedron_rules.cpp
				grid_objects/rule_util.cpp
				multi_grid.cpp
				lib_grid_messages.cpp
				)


if(PARALLEL)
	# add parallelization
	set(srcParallelization	parallelization/broadcast.cpp
							parallelization/distribution.cpp
							parallelization/distributed_grid.cpp
							parallelization/gather_grid.cpp
							parallelization/parallelization_util.cpp
							parallelization/parallel_grid_layout.cpp
							parallelization/load_balancer.cpp
							parallelization/load_balancer_util.cpp
							parallelization/deprecated/load_balancing.cpp
							parallelization/partitioner_dynamic_bisection.cpp
							parallelization/parallel_refinement/parallel_global_fractured_media_refiner.cpp
							parallelization/parallel_refinement/parallel_global_subdivision_refiner.cpp
							parallelization/parallel_refinement/parallel_hanging_node_refiner_multi_grid.cpp
							parallelization/parallel_refinement/parallel_hnode_adjuster.cpp)
	
else(PARALLEL)
	set(srcParallelization	parallelization/deprecated/load_balancing.cpp
							parallelization/parallel_grid_layout.cpp)
endif(PARALLEL)

################################################
# EXTERNALS - liblgm
set(srcLgm
	file_io/externals/src/lgm/src/lgm.c
	file_io/externals/src/lgm/src/lgm_error.c
	file_io/externals/src/lgm/src/lgm_info.c
	file_io/externals/src/lgm/src/lgm_parser.c
	file_io/externals/src/lgm/src/lgm_writer.c
)

################################################
# EXTERNALS - libng
set(srcNg
	file_io/externals/src/ng/src/ng.c
	file_io/externals/src/ng/src/ng_error.c
	file_io/externals/src/ng/src/ng_info.c
	file_io/externals/src/ng/src/ng_parser.c
	file_io/externals/src/ng/src/ng_writer.c
)

set(srcTokstream
	file_io/externals/src/tokstream/tokstream.c
)

#the grid library					

set(Sources         ${srcLibGrid}
                    ${srcParallelization}
					${srcGrid} 
					${srcAlgorithms} ${srcFileIO}
					${srcTokstream} ${srcTetgen} ${srcLgm} ${srcNg})

if(BUILD_ONE_LIB)
	EXPORTSOURCES(ugbase/lib_grid ${Sources})
else(BUILD_ONE_LIB)
	add_library(grid STATIC	${Sources})
endif(BUILD_ONE_LIB)

#for debug only
#add_library(grid ${srcGrid} ${srcLibGrid})

################################################
# created by Sebastian Reiter
# s.b.reiter@googlemail.com
#
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
				grid/neighborhood.cpp)
				
set(srcAlgorithms	algorithms/debug_util.cpp
					algorithms/subset_util.cpp
					algorithms/selection_util.cpp
					algorithms/serialization.cpp
					algorithms/grid_statistics.cpp
					algorithms/polychain_util.cpp
					algorithms/problem_detection_util.cpp
					algorithms/remeshing/edge_length_adjustment.cpp
					algorithms/remeshing/simple_grid.cpp
					algorithms/remeshing/grid_adaption.cpp
					algorithms/geom_obj_util/edge_util.cpp
					algorithms/geom_obj_util/face_util.cpp
					algorithms/geom_obj_util/misc_util.cpp
					algorithms/geom_obj_util/vertex_util.cpp
					algorithms/geom_obj_util/volume_util.cpp
					algorithms/grid_generation/icosahedron.cpp
					algorithms/grid_generation/tetrahedralization.cpp
					algorithms/grid_generation/triangle_fill.cpp
					algorithms/grid_generation/triangle_fill_sweep_line.cpp
					algorithms/extrusion/extrude.cpp
					algorithms/extrusion/cylinder_extrusion.cpp
					algorithms/extrusion/expand_layers.cpp
					algorithms/refinement/ref_mark_adjusters/std_hnode_adjuster.cpp
					algorithms/refinement/ref_mark_adjusters/mg_hnode_adjuster.cpp
					algorithms/refinement/refinement_projectors/misc_refinement_projectors.cpp
					algorithms/refinement/refinement_projectors/fractal_projector.cpp
					algorithms/refinement/adaptive_regular_mg_refiner.cpp
					algorithms/refinement/fractured_media_refiner.cpp
					algorithms/refinement/hanging_node_refiner_base.cpp
					algorithms/refinement/hanging_node_refiner_grid.cpp
					algorithms/refinement/hanging_node_refiner_multi_grid.cpp
					algorithms/refinement/refiner_interface.cpp
					# (unmaintained) algorithms/refinement/hanging_node_refiner_2d_irn.cpp
					# (unmaintained) algorithms/refinement/multi_grid_refiner.cpp
					algorithms/refinement/regular_refinement.cpp
					algorithms/refinement/global_multi_grid_refiner.cpp
					algorithms/refinement/global_fractured_media_refiner.cpp
					algorithms/subdivision/subdivision_loop.cpp
					algorithms/subdivision/subdivision_rules_piecewise_loop.cpp
					algorithms/unit_tests/check_associated_elements.cpp)
					
set(srcFileIO	file_io/file_io_art.cpp
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
	set(srcParallelization	parallelization/distribution.cpp
							parallelization/distributed_grid.cpp
							parallelization/parallelization_util.cpp
							parallelization/parallel_grid_layout.cpp
							parallelization/load_balancer.cpp
							parallelization/load_balancer_util.cpp
							parallelization/load_balancing.cpp
							parallelization/partitioner_dynamic_bisection.cpp
							parallelization/parallel_refinement/parallel_global_fractured_media_refiner.cpp
							parallelization/parallel_refinement/parallel_hanging_node_refiner_multi_grid.cpp
							parallelization/parallel_refinement/parallel_hnode_adjuster.cpp)
	if(PARMETIS)
		set(srcParallelization	${srcParallelization}
								parallelization/partitioner_parmetis.cpp)
	endif(PARMETIS)
	
else(PARALLEL)
	set(srcParallelization	parallelization/load_balancing.cpp
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

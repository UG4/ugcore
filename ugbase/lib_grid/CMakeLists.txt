################################################
# created by Sebastian Reiter
# s.b.reiter@googlemail.com
#
# lib_grid
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_GRID)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcGrid		grid/grid.cpp
				grid/geometric_base_objects.cpp
				grid/grid_connection_managment.cpp
				grid/geometric_object_collection.cpp
				grid/grid_util.cpp)
				
set(srcAlgorithms	algorithms/subset_util.cpp
					algorithms/selection_util.cpp
					algorithms/serialization.cpp
					algorithms/grid_statistics.cpp
					algorithms/extruder_util.cpp
					algorithms/remeshing/edge_length_adjustment.cpp
					algorithms/remeshing/simple_grid.cpp
					algorithms/remeshing/grid_adaption.cpp
					algorithms/graph/neighbourhood.cpp
					algorithms/geom_obj_util/edge_util.cpp
					algorithms/geom_obj_util/face_util.cpp
					algorithms/geom_obj_util/misc_util.cpp
					algorithms/geom_obj_util/vertex_util.cpp
					algorithms/geom_obj_util/volume_util.cpp
					algorithms/grid_generation/tetrahedralization.cpp
					algorithms/extrusion/extrusion.cpp
					algorithms/refinement/hanging_node_refiner.cpp
					algorithms/refinement/multi_grid_refiner.cpp
					algorithms/refinement/regular_refinement.cpp
					algorithms/refinement/global_multi_grid_refiner.cpp
					algorithms/subdivision/subdivision_loop.cpp)
					
set(srcFileIO	file_io/file_io_art.cpp
				file_io/file_io_dump.cpp
				file_io/file_io_lgb.cpp
				file_io/file_io_lgm.cpp
				file_io/file_io_ng.cpp
				file_io/file_io_obj.cpp
				file_io/file_io_tetgen.cpp
				file_io/file_io_txt.cpp
				file_io/file_io_ug.cpp
				file_io/file_io.cpp)
								
set(srcLibGrid	tools/subset_handler_interface.cpp
				tools/subset_handler_grid.cpp
				tools/subset_handler_multi_grid.cpp
				tools/selector_interface.cpp
				tools/selector_grid.cpp
				tools/selector_multi_grid.cpp
				tools/selector_grid_elem.cpp
				tools/marker_points.cpp
				geometric_objects.cpp
				common_attachments.cpp
				multi_grid.cpp
				)


if(COMPILE_WITH_PCL)
	# add parallelization
	set(srcParallelization	algorithms/parallelization/distribution_util.cpp
							algorithms/parallelization/grid_distribution.cpp
							algorithms/parallelization/distributed_grid.cpp
							algorithms/parallelization/parallel_refinement/parallel_multi_grid_refiner.cpp
							algorithms/parallelization/parallel_refinement/parallel_global_multi_grid_refiner.cpp)

else(COMPILE_WITH_PCL)
	# give a warning that pcl will not be compiled.
	set(srcParallelization )
endif(COMPILE_WITH_PCL)

################################################
# EXTERNALS - TetGen
add_definitions(-DTETLIBRARY)

set_source_files_properties(externals/src/tetgen/predicates.cxx
	PROPERTIES COMPILE_FLAGS -O0)

set_source_files_properties(externals/src/tetgen/tetgen.cxx
	PROPERTIES COMPILE_FLAGS -w)

set(srcTetgen
	externals/src/tetgen/tetgen.cxx
	externals/src/tetgen/predicates.cxx
)

################################################
# EXTERNALS - liblgm
set(srcLgm
	externals/src/lgm/src/lgm.c
	externals/src/lgm/src/lgm_error.c
	externals/src/lgm/src/lgm_info.c
	externals/src/lgm/src/lgm_parser.c
	externals/src/lgm/src/lgm_writer.c
	externals/src/lgm/src/tokstream/tokstream.c
)

################################################
# EXTERNALS - libng
set(srcNg
	externals/src/ng/src/ng.c
	externals/src/ng/src/ng_error.c
	externals/src/ng/src/ng_info.c
	externals/src/ng/src/ng_parser.c
	externals/src/ng/src/ng_writer.c
	externals/src/ng/src/tokstream/tokstream.c
)

#the real library					
add_library(grid 	${srcGrid} ${srcLibGrid}
					${srcParallelization}
					${srcAlgorithms} ${srcFileIO}
					${srcTetgen} ${srcLgm} ${srcNg})

#for debug only
#add_library(grid ${srcGrid} ${srcLibGrid})

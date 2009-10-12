################################################
# created by Sebastian Reiter
# s.b.reiter@googlemail.com
#
# lib_grid
################################################

cmake_minimum_required(VERSION 2.6)

project(P_LIB_GRID)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcGrid		grid/grid_connection_managment.cpp
				grid/geometric_object_collection.cpp
				grid/grid_util.cpp
				grid/grid.cpp)
				
set(srcAdvanced	advanced/extrusion/extrusion.cpp
				advanced/refinement/hanging_node_refiner.cpp
				advanced/refinement/multi_grid_refiner.cpp)

set(srcAlgorithms	algorithms/subset_util.cpp
					algorithms/geom_obj_util/edge_util.cpp
					algorithms/geom_obj_util/face_util.cpp
					algorithms/geom_obj_util/misc_util.cpp
					algorithms/geom_obj_util/vertex_util.cpp
					algorithms/grid_generation/tetrahedralization.cpp)
					
set(srcFileIO	file_io/file_io_dump.cpp
				file_io/file_io_lgm.cpp
				file_io/file_io_ng.cpp
				file_io/file_io_obj.cpp
				file_io/file_io_tetgen.cpp
				file_io/file_io_txt.cpp
				file_io/file_io_ug.cpp
				file_io/file_io.cpp)
								
set(srcLibGrid	common_attachments.cpp
				geometric_objects.cpp
				multi_grid.cpp
				selector.cpp
				subset_handler.cpp)

				
################################################
# EXTERNALS - TetGen
add_definitions(-DTETLIBRARY)

set_source_files_properties(externals/src/tetgen/predicates.cxx
	PROPERTIES COMPILE_FLAGS -O0)

set_source_files_properties(externals/src/tetgen/tetgen.cxx
	PROPERTIES COMPILE_FLAGS -Wno-write-strings)

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
					
add_library(grid 	${srcGrid} ${srcAdvanced} ${srcAlgorithms}
					${srcFileIO} ${srcLibGrid}
					${srcTetgen} ${srcLgm} ${srcNg})


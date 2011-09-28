cmake_minimum_required(VERSION 2.6)

################################################
# node-tree
################################################

project(P_NODE_TREE)

include_directories(${P_UGBASE_SOURCE_DIR})

set(srcNodeTree	boxed_group_node.cpp
				collision_edges_node.cpp
				collision_triangles_node.cpp
				collision_tree_root_node.cpp
				group_node.cpp
				traverser.cpp
				octree.cpp)

set(srcTraversers	traversers/traverser_collision_tree.cpp
					traversers/traverser_project_point.cpp
					traversers/traverser_intersect_faces.cpp)



set(Sources     ${srcNodeTree} ${srcTraversers})


if(BUILD_ONE_LIB)
	ExportSources("common/node_tree" ${Sources})
else(BUILD_ONE_LIB)
	add_library(node_tree STATIC ${Sources})
endif(BUILD_ONE_LIB)

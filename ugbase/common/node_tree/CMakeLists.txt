# Copyright (c) 2010-2012:  G-CSC, Goethe University Frankfurt
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
	ExportSources("ugbase/common/node_tree" ${Sources})
else(BUILD_ONE_LIB)
	add_library(node_tree STATIC ${Sources})
endif(BUILD_ONE_LIB)

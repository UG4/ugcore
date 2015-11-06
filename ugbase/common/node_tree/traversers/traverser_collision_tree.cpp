/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include "traverser_collision_tree.h"
#include "../node_tree.h"

namespace ug{
namespace node_tree
{

Traverser_CollisionTree::Traverser_CollisionTree()
{
//	register handler functions
	register_handler_function(OC_COLLISION_TREE_ROOT_NODE,
							&Traverser_CollisionTree::handle_collision_tree_root);
	register_handler_function(OC_COLLISION_EDGES_NODE,
							&Traverser_CollisionTree::handle_collision_edges);
	register_handler_function(OC_COLLISION_TRIANGLES_NODE,
							&Traverser_CollisionTree::handle_collision_triangles);
}

Traverser_CollisionTree::~Traverser_CollisionTree()
{
}

void Traverser_CollisionTree::
handle_collision_tree_root(CollisionTreeRootNode* colTreeRootNode)
{
	//std::cout << "handling collision tree root" << std::endl;

//	put the rootNode on top of the stack.
//	All CollisionEdgeNodes children will index it's points
	m_stackRootNodes.push(colTreeRootNode);

//	traverse the node
	handle_boxed_group(colTreeRootNode);

//	pop the rootNode from the stack
	m_stackRootNodes.pop();
}

void Traverser_CollisionTree::
handle_collision_edges(CollisionEdgesNode* colEdgesNode)
{
	//std::cout << "handling collision edges" << std::endl;
}

void Traverser_CollisionTree::
handle_collision_triangles(CollisionTrianglesNode* colTrisNode)
{
	//std::cout << "handling collision edges" << std::endl;
}

CollisionTreeRootNode* Traverser_CollisionTree::
get_current_root_node()
{
	assert(!m_stackRootNodes.empty() && "Make sure that your collision tree"
			" contains a CollisionTreeRootNode!");
	return m_stackRootNodes.top();
}

}//	end of namespace node_tree
}//	end of namespace ug

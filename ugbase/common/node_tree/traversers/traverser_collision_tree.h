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

#ifndef __H__UG__NODE_TREE__TRAVERSER_COLLISION_TREE__
#define __H__UG__NODE_TREE__TRAVERSER_COLLISION_TREE__

#include <stack>
#include "../traverser.h"
#include "../collision_element_info.h"

namespace ug{
namespace node_tree
{

class CollisionTreeRootNode;
class CollisionEdgesNode;
class CollisionTrianglesNode;

////////////////////////////////////////////////////////////////////////
//	Traverser_CollisionTree
///	Enhances the Traverser base-class by methods to traverse a collision tree.
/**
 * Derive from this class if you want to create a traverser that traverses
 * a collision-tree.
 *
 * This class does nothing more than registering its callback-methods at the
 * base-traverser.
 *
 * The methods here do not do much. handle_collision_tree_root calls its
 * handle_boxed_group on its node. The other methods are empty.
 */
class Traverser_CollisionTree : public Traverser
{
	public:
		Traverser_CollisionTree();
		virtual ~Traverser_CollisionTree();

	protected:
		virtual void handle_collision_tree_root(CollisionTreeRootNode* colTreeRootNode);
		virtual void handle_collision_edges(CollisionEdgesNode* colEdgesNode);
		virtual void handle_collision_triangles(CollisionTrianglesNode* colTrisNode);

		CollisionTreeRootNode* get_current_root_node();

	private:
		std::stack<CollisionTreeRootNode*>	m_stackRootNodes;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif

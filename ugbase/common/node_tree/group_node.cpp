/*
 * Copyright (c) 2007-2015:  Sebastian Reiter
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

#include "group_node.h"

namespace ug{
namespace node_tree
{

SPGroupNode GroupNode::create()
{
	GroupNode* node = new GroupNode;
	node->m_objectCode = OC_GROUP_NODE;
	return SPGroupNode(node);
}

GroupNode::GroupNode()
{
}

GroupNode::~GroupNode()
{
}

void GroupNode::clear()
{
	vChildren.clear();
}

void GroupNode::add_child(SPNode node)
{
	vChildren.push_back(node);
}

void GroupNode::remove_child(SPNode node)
{
	for(std::vector<SPNode>::iterator iter = vChildren.begin();
		iter != vChildren.end(); iter++)
	{
		SPNode& spNode = *iter;
		if(spNode == node)
		{
			vChildren.erase(iter);
			break;
		}
	}
}

int GroupNode::num_children()
{
	return vChildren.size();
}

SPNode GroupNode::get_child(int index)
{
	return vChildren[index];
}

}//	end of namespace node_tree
}//	end of namespace ug

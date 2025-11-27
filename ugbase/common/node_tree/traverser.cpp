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

#include "traverser.h"
#include "node_tree.h"
#include <iostream>

namespace ug{
namespace node_tree
{

Traverser::Traverser()
{
//	initialise the handler-function array to reduce the number of resizes
	m_vHandlerFuncs.resize(OC_NODES_END, 0);

//	register handler functions
	register_handler_function(OC_GROUP_NODE, &Traverser::handle_group);
	register_handler_function(OC_BOXED_GROUP_NODE, &Traverser::handle_boxed_group);
}


bool Traverser::handler_function_registered(unsigned int oc) const {
	// the first expression has been always true for unsigned int (AV)
	if(/*oc >= 0 &&*/ oc < m_vHandlerFuncs.size())
	{
		if(m_vHandlerFuncs[oc] != 0)
			return true;
	}
	return false;
}

void Traverser::apply(SPNode& node)
{
	traverse_object(node.get());
}

void Traverser::traverse_object(Object* obj)
{
	unsigned int oc = obj->getObjectCode();
	if(handler_function_registered(oc))
	{
		(this->*(m_vHandlerFuncs[oc]))(obj);
	}
}

void Traverser::handle_group(GroupNode* group)
{
	//std::cout << "handling group" << std::endl;

//	traverse all children of the group
	int numChildren = group->num_children();
	for(int i = 0; i < numChildren; ++i)
		traverse_object(group->get_child(i).get());
}

void Traverser::handle_boxed_group(BoxedGroupNode* boxedGroup)
{
	//std::cout << "handling boxed group" << std::endl;

//	traverse the group
	handle_group(boxedGroup);
}

}//	end of namespace node_tree
}//	end of namespace ug

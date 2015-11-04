//	traverser.cpp
//	created by Sebastian Reiter y07 m12 d4
//	s.b.reiter@googlemail.com

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

Traverser::~Traverser()
{
}

bool Traverser::handler_function_registered(unsigned int oc)
{
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

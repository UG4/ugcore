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

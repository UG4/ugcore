#include "boxed_group_node.h"

namespace ug{
namespace node_tree
{

SPBoxedGroupNode BoxedGroupNode::create()
{
	BoxedGroupNode* node = new BoxedGroupNode;
	node->m_objectCode = OC_BOXED_GROUP_NODE;
	return SPBoxedGroupNode(node);
}

BoxedGroupNode::BoxedGroupNode()
{
}

BoxedGroupNode::~BoxedGroupNode()
{
}

void BoxedGroupNode::set_box(const vector3& minCorner,
							 const vector3& maxCorner)
{
	m_vMin = minCorner;
	m_vMax = maxCorner;
}

const vector3&	BoxedGroupNode::min_corner() const
{
	return m_vMin;
}

const vector3&	BoxedGroupNode::max_corner() const
{
	return m_vMax;
}

}//	end of namespace node_tree
}//	end of namespace ug

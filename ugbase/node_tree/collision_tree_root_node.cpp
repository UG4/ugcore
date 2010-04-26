//	collision_tree_root_node.cpp
//	created by Sebastian Reiter y07 m12 d5
//	s.b.reiter@googlemail.com

#include "collision_tree_root_node.h"

namespace ug{
namespace node_tree
{

////////////////////////////////////////////////////////////////////////////////////////////////
SPCollisionTreeRootNode CollisionTreeRootNode::create()
{
	CollisionTreeRootNode* node = new CollisionTreeRootNode;
	node->m_objectCode = OC_COLLISION_TREE_ROOT_NODE;
	return SPCollisionTreeRootNode(node);
}

////////////////////////////////////////////////////////////////////////////////////////////////
CollisionTreeRootNode::CollisionTreeRootNode()
{
}

CollisionTreeRootNode::~CollisionTreeRootNode()
{
}

void CollisionTreeRootNode::clear_points()
{
	m_vPoints.clear();
}

int CollisionTreeRootNode::num_points()
{
	return m_vPoints.size();
}

void CollisionTreeRootNode::add_points(vector3* pPoints, int numPoints)
{
	for(int i = 0; i < numPoints; ++i)
		m_vPoints.push_back(pPoints[i]);
}

const vector3& CollisionTreeRootNode::get_point(int index) const
{
	return m_vPoints[index];
}

const vector3* CollisionTreeRootNode::get_points() const
{
	return &m_vPoints.front();
}

}//	end of namespace node_tree
}//	end of namespace ug

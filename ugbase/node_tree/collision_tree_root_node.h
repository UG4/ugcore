//	collision_tree_root_node.h
//	created by Sebastian Reiter y07 m12 d4
//	s.b.reiter@googlemail.com

#ifndef __H__UG__NODE_TREE__COLLISION_TREE_ROOT_NODE__
#define __H__UG__NODE_TREE__COLLISION_TREE_ROOT_NODE__

#include "boxed_group_node.h"
#include <vector>

namespace ug{
namespace node_tree
{

class CollisionTreeRootNode;
////////////////////////////////////////////////////////////////////////
///	the smartpointer which will be used to encapsulate the node
typedef SmartPtr<CollisionTreeRootNode> SPCollisionTreeRootNode;

////////////////////////////////////////////////////////////////////////
//	CollisionTreeRootNode
///	A group node featuring a bounding box and a set of points.
/**
* Points are indexed by the subsidiary leaf elements.
*/
class CollisionTreeRootNode : public BoxedGroupNode
{
	public:
		static SPCollisionTreeRootNode create();

		virtual ~CollisionTreeRootNode();

		virtual void clear_points();
		virtual int num_points();
		virtual void add_points(vector3* pPoints, int numPoints);
		virtual const vector3& get_point(int index) const;
		virtual const vector3* get_points() const;

	protected:
		CollisionTreeRootNode();

	protected:
		typedef std::vector<vector3> PointVec;
		PointVec	m_vPoints;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif

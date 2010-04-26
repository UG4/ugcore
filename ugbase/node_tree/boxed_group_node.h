//	boxed_group_node.h
//	created by Sebastian Reiter y07 m12 d4
//	s.b.reiter@googlemail.com

#ifndef __H__UG__NODE_TREE__BOXED_GROUP_NODE__
#define __H__UG__NODE_TREE__BOXED_GROUP_NODE__

#include "group_node.h"
#include "common/math/ugmath.h"

namespace ug{
namespace node_tree
{

class BoxedGroupNode;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
typedef SmartPtr<BoxedGroupNode> SPBoxedGroupNode;

////////////////////////////////////////////////////////////////////////
//	BoxedGroupNode
///	A group node featuring a bounding box
/**
...
*/
class BoxedGroupNode : public GroupNode
{
	public:
		static SPBoxedGroupNode create();

		virtual ~BoxedGroupNode();

		virtual void set_box(const vector3& minCorner,
							 const vector3& maxCorner);

		virtual const vector3& min_corner() const;
		virtual const vector3& max_corner() const;

	protected:
		BoxedGroupNode();

	protected:
		vector3 m_vMin;
		vector3 m_vMax;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif

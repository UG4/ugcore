//	group_node.h
//	created by Sebastian Reiter 15. November 2007
//	s.b.reiter@googlemail.com

#ifndef __H__UG__NODE_TREE__GROUP_NODE__
#define __H__UG__NODE_TREE__GROUP_NODE__

#include <vector>
#include "node.h"

namespace ug{
namespace node_tree
{

class GroupNode;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
typedef SmartPtr<GroupNode> SPGroupNode;

////////////////////////////////////////////////////////////////////////
//	GroupNode
///	You can group nodes using the GroupNode
/**
...
*/
class GroupNode : public Node
{
	public:
		static SPGroupNode create();

		virtual ~GroupNode();
		
		virtual void clear();
		virtual void add_child(SPNode node);
		virtual void remove_child(SPNode node);
		
		virtual int num_children();
		virtual SPNode get_child(int index);

	protected:
		GroupNode();

	protected:
		std::vector<SPNode>	vChildren;
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif

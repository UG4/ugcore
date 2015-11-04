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

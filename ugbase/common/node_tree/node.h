#ifndef __H__UG__NODE_TREE__NODE__
#define __H__UG__NODE_TREE__NODE__

#include "object.h"
#include "common/util/smart_pointer.h"


////////////////////////////////////////////////////////////////////////
//	SPNODE
///	casts any type to SPNode.
/**
 * This macro should only be applied to SmartPointers,
 * which encapsulate a derivative of Node.
 */
//not needed any more, if the smart-pointer template constructor works
//TODO: remove this makro.
#define SPNODE(node) *((ng::SPNode*)&node)


namespace ug{
namespace node_tree
{

class Node;

////////////////////////////////////////////////////////////////////////
///	the smartpointer used to encapsulate the node
typedef SmartPtr<Node> SPNode;

////////////////////////////////////////////////////////////////////////
//	Node
///	Nodes serve as base-objects for items of wich a scene-graph consists.
/**
 * Nodes can not be directly instanced. Derivatives of Node should
 * feature a static create() function, through wich the user can recieve
 * an instance (encapsulated by the smart pointer SPNode).
 * This concept should be shared by all derived nodes.
 */
class Node : public Object
{
	friend class Traverser;
	
	public:
		virtual ~Node()	{}
		
	protected:
		Node()	{}
};

}//	end of namespace node_tree
}//	end of namespace ug

#endif

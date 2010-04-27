//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d27

#include "traverser_collision_tree.h"
#include "../node_tree.h"

namespace ug{
namespace node_tree
{

Traverser_CollisionTree::Traverser_CollisionTree()
{
//	register handler functions
	register_handler_function(OC_COLLISION_TREE_ROOT_NODE,
							&Traverser_CollisionTree::handle_collision_tree_root);
	register_handler_function(OC_COLLISION_EDGES_NODE,
							&Traverser_CollisionTree::handle_collision_edges);
	register_handler_function(OC_COLLISION_TRIANGLES_NODE,
							&Traverser_CollisionTree::handle_collision_triangles);
}

Traverser_CollisionTree::~Traverser_CollisionTree()
{
}

void Traverser_CollisionTree::
handle_collision_tree_root(CollisionTreeRootNode* colTreeRootNode)
{
	//std::cout << "handling collision tree root" << std::endl;

//	handle the boxed group
	handle_boxed_group(colTreeRootNode);
}

void Traverser_CollisionTree::
handle_collision_edges(CollisionEdgesNode* colEdgesNode)
{
	//std::cout << "handling collision edges" << std::endl;
}

void Traverser_CollisionTree::
handle_collision_triangles(CollisionTrianglesNode* colTrisNode)
{
	//std::cout << "handling collision edges" << std::endl;
}

}//	end of namespace node_tree
}//	end of namespace ug

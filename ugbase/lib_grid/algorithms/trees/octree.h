//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d26

#ifndef __H__LIB_GRID__OCTREE__
#define __H__LIB_GRID__OCTREE__

#include "node_tree/node_tree.h"
#include "lib_grid/lg_base.h"

namespace ug
{



////////////////////////////////////////////////////////////////////////
///	Creates an Octree from a list of edges or triangles.
/**
 * The elements to which the iterators point have to be derivates of
 * ug::EdgeBase or ug::Face and have to represent edges or triangles.
 * Furthermore all elements between elemsBegin and elemsBegin have to 
 * be members of the specified grid.
 *
 * Each element in the resulting octree will store a pointer to the
 * associated edge or triangle of the grid. You may access this pointer
 * using the CollisionEdgesNode::get_edge_id or
 * CollisionTrianglesNode::get_triangle_id method to retreive it.
 *
 * This method internally uses the ug::node_tree::CreateOctree method.
 */
template <class TIterator>
node_tree::SPCollisionTreeRootNode
CreateOctree(Grid& grid, TIterator elemsBegin, TIterator elemsEnd,
			int maxDepth, int elemThreshold, bool bLoose,
			APosition& aPos = aPosition);

							   
}//	end of namespace

////////////////////////////////
//	include implementation
#include "octree_impl.hpp"

#endif

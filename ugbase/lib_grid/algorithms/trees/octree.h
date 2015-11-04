#ifndef __H__LIB_GRID__OCTREE__
#define __H__LIB_GRID__OCTREE__

#include "common/node_tree/node_tree.h"
#include "lib_grid/lg_base.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_trees
///	@{

////////////////////////////////////////////////////////////////////////
//	SPOctree
/**	This typedef is only intended to increase comfort for callers.
 *	The Octree can be passed to all traversers to which the original
 *	CollisionTreeRootNode could be passed too.
 *	\sa node_tree::SPCollisionTreeRootNode
 */
typedef node_tree::SPCollisionTreeRootNode SPOctree;

////////////////////////////////////////////////////////////////////////
///	Creates an Octree from a list of edges or triangles.
/**
 * This method returns a SmartPointer to an Octree. You may access
 * all public members using the -> operator.
 *
 * The elements to which the iterators point have to be derivates of
 * ug::Edge or ug::Face and have to represent edges or triangles.
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
SPOctree
CreateOctree(Grid& grid, TIterator elemsBegin, TIterator elemsEnd,
			int maxDepth, int elemThreshold, bool bLoose,
			APosition& aPos = aPosition);

/// @}
							   
}//	end of namespace

////////////////////////////////
//	include implementation
#include "octree_impl.hpp"

#endif

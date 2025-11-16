/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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
/**	This type definition is only intended to increase comfort for callers.
 *	The Octree can be passed to all traversers to which the original
 *	CollisionTreeRootNode could be passed too.
 *	\sa node_tree::SPCollisionTreeRootNode
 */
using SPOctree = node_tree::SPCollisionTreeRootNode;

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

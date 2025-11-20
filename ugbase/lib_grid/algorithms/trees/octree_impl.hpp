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

#ifndef __H__LIB_GRID__OCTREE_IMPL__
#define __H__LIB_GRID__OCTREE_IMPL__

#include <vector>
#include "lib_grid/lg_base.h"
#include "common/node_tree/octree.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CreateOctree
template <typename TIterator>
SPOctree
CreateOctree(Grid& grid, TIterator elemsBegin, TIterator elemsEnd,
			int maxDepth, int elemThreshold, bool bLoose,
			APosition& aPos)
{
	const char* logMsgPrefix = "  CreateOctree: ";

	if(elemsBegin == elemsEnd)
		return node_tree::SPCollisionTreeRootNode(nullptr);

//	access the position attachment of the grid.	
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG(logMsgPrefix << "vertex-attachment missing: aPos\n");
		return node_tree::SPCollisionTreeRootNode(nullptr);
	}
	
	Grid::VertexAttachmentAccessor aaPos(grid, aPos);
	
//	we have to attach an int-attachment to the vertices of the grid.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor aaInt(grid, aInt);
	
//	all elements have to have the same amount of vertices
	const size_t numVertices = (*elemsBegin)->num_vertices();
	
//	initialise all entries of aInt to -1
//	since we want to allow to construct an OctTree even for small subsets
//	of a grid, we only want to touch vertices that are important.
	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		for(size_t i = 0; i < numVertices; ++i)
			aaInt[(*iter)->vertex(i)] = -1;
	}
	
//	create the vElems and the vPoints arrays.
	std::vector<vector3> vPoints;
	std::vector<int> vElems;
	std::vector<node_tree::CollisionElementID> vElemIDs;

/* If you want to match the indices that are referenced by the oct-trees
 * triangles with the order of the grids vertices, than you can uncomment
 * the following code. Please note, that this code is slower than the
 * original code, if only a tiny subset of the grid is sorted into the
 * tree.
int tmpCounter = 0;
for(VertexIterator iter = grid.vertices_begin(); iter != grid.vertices_end(); ++iter){
	aaInt[*iter] = tmpCounter++;
	vPoints.push_back(aaPos[*iter]);
}
*/

	for(TIterator iter = elemsBegin; iter != elemsEnd; ++iter){
		typename TIterator::value_type elem = *iter;
		
	//	insert the indices into vElems and build vPoints at the same time
		for(size_t i = 0; i < numVertices; ++i){
			Vertex* v = elem->vertex(i);
			if(aaInt[v] == -1){
				aaInt[v] = (int)vPoints.size();
				vPoints.push_back(aaPos[v]);
			}
			
			vElems.push_back(aaInt[v]);
		}
		
	//	insert the element-id (the pointer to the associated grid element)
		vElemIDs.push_back(elem);
	}
	
//	arrays have been created. remove obsolete attachments
	grid.detach_from_vertices(aInt);
	
//	call the octree creation algorithm
	return node_tree::CreateOctree(&vPoints.front(), vPoints.size(),
							&vElems.front(), vElems.size(), numVertices,
							&vElemIDs.front(), maxDepth, elemThreshold, bLoose);
}

}//	end of namespace

#endif

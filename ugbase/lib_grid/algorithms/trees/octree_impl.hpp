//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m04 d26

#ifndef __H__LIB_GRID__OCTREE_IMPL__
#define __H__LIB_GRID__OCTREE_IMPL__

#include <vector>
#include "lib_grid/lg_base.h"
#include "common/node_tree/octree.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	CreateOctree
template <class TIterator>
SPOctree
CreateOctree(Grid& grid, TIterator elemsBegin, TIterator elemsEnd,
			int maxDepth, int elemThreshold, bool bLoose,
			APosition& aPos)
{
	const char* logMsgPrefix = "  CreateOctree: ";

	if(elemsBegin == elemsEnd)
		return node_tree::SPCollisionTreeRootNode(NULL);

//	access the position attachment of the grid.	
	if(!grid.has_vertex_attachment(aPos)){
		UG_LOG(logMsgPrefix << "vertex-attachment missing: aPos\n");
		return node_tree::SPCollisionTreeRootNode(NULL);
	}
	
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	
//	we have to attach an int-attachment to the vertices of the grid.
	AInt aInt;
	grid.attach_to_vertices(aInt);
	Grid::VertexAttachmentAccessor<AInt> aaInt(grid, aInt);
	
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
	std::vector<vector3> 	vPoints;
	std::vector<int>		vElems;
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

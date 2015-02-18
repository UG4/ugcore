// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_mark_util_impl
#define __H__UG_mark_util_impl

#include "lib_grid/algorithms/geom_obj_util/edge_util.h"
namespace ug{

template <class TEdgeIterator>
void MarkCreaseEdges(Grid& grid, ISubsetHandler& sh,
					TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
					int subsetIndex, number angle,
					APosition& aPos,
					ANormal* paFaceNormal)
{

//	get the position accessor
	if(!grid.has_vertex_attachment(aPos))
		return;
	
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);
	
//	we'll store face normals in those vectors:
	vector3 n[2];
	
//	associated faces are stored in this array
	Face* f[2];
	
//	all dot-products between normals lower than minDot mark a crease.
	number minDot = cos(angle * 3.14159265 / 180.f);
	
//	iterate through the edges
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
		Edge* e = *iter;
	//	get the associated faces
	//	all edges that do not have exactly 2 associated edges
	//	are regarded as seam-edges
		if(GetAssociatedFaces(f, grid, e, 2) == 2){
		//	get the normals of the associated faces
			CalculateNormal(n[0], f[0], aaPos);
			CalculateNormal(n[1], f[1], aaPos);
		//	if the dot-product is lower than minDot, then the edge is a crease edge.
			if(VecDot(n[0], n[1]) < minDot)
				sh.assign_subset(e, subsetIndex);
		}
		else{
			sh.assign_subset(e, subsetIndex);
		}
	}
}

template <class TVertexIterator, class TAPosition>
void MarkCorners(Grid& grid, ISubsetHandler& sh,
				 TVertexIterator vrtsBegin, TVertexIterator vrtsEnd,
				 Grid::edge_traits::callback cbPathEdge,
				 int subsetIndex, number angle,
				 TAPosition& aPos)
{
	typedef typename TAPosition::ValueType	vector_t;

	if(!grid.has_vertex_attachment(aPos))
		return;
	
	Grid::VertexAttachmentAccessor<APosition> aaPos(grid, aPos);

	Grid::edge_traits::secure_container	edges;
	vector_t	dir[2];
	Edge*		pathEdge[2];

	number dotThreshold = cos(deg_to_rad(angle));

	for(TVertexIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		Vertex* vrt = *iter;
		grid.associated_elements(edges, vrt);
		int numPathEdges = 0;

		for(size_t i = 0; i < edges.size(); ++i){
			Edge* e = edges[i];
			if(cbPathEdge(e)){
				if(numPathEdges < 2)
					pathEdge[numPathEdges] = e;
				++numPathEdges;
			}
		}

		if(numPathEdges == 2){
			for(size_t i = 0; i < 2; ++i){
				Vertex* cv = GetConnectedVertex(pathEdge[i], vrt);
				VecSubtract(dir[i], aaPos[cv], aaPos[vrt]);
				VecNormalize(dir[i], dir[i]);
			}

			dir[1] *= -1;
			if(VecDot(dir[0], dir[1]) < dotThreshold)
				sh.assign_subset(vrt, subsetIndex);
		}
		else if(numPathEdges > 0){
			sh.assign_subset(vrt, subsetIndex);
		}
	}
}

}//	end of namespace

#endif	//__H__mark_util_impl

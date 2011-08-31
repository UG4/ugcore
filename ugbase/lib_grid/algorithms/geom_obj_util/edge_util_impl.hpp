// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y08 m12 d16

#ifndef __H__LIB_GRID__EDGE_UTIL_IMPL__
#define __H__LIB_GRID__EDGE_UTIL_IMPL__

//#include "edge_util.h"
#include <stack>
#include "lib_grid/grid/grid_util.h"
#include "vertex_util.h"
namespace ug
{

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
inline number EdgeLengthSq(EdgeBase* e, TAAPosVRT& aaPos)
{
	return VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
inline number EdgeLength(EdgeBase* e, TAAPosVRT& aaPos)
{
	return VecDistance(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
}

////////////////////////////////////////////////////////////////////////
//	SplitEdge
//	see edge_operations.h for detailed description
template<class TVertex>
TVertex* SplitEdge(Grid& grid, EdgeBase* e, bool bConservative)
{
	return SplitEdge<TVertex>(grid, grid, e, NULL, bConservative);
}

////////////////////////////////////////////////////////////////////////
//	SplitEdge
//	see edge_operations.h for detailed description
template<class TVertex>
TVertex* SplitEdge(Grid& destGrid, Grid& srcGrid, EdgeBase* e,
						AVertexBase* paAssociatedVertices,
						bool bConservative)
{
	TVertex* newVertex;
	if(&destGrid == &srcGrid)
		newVertex = *destGrid.create<TVertex>(e);
	else
		newVertex = *destGrid.create<TVertex>();

	if(CreateEdgeSplitGeometry(destGrid, srcGrid, e, newVertex, paAssociatedVertices))
	{
		if(!bConservative)
		{
		//	erase unused elements.
			if(!srcGrid.option_is_enabled(FACEOPT_AUTOGENERATE_EDGES))
			{
			//	we have to erase the faces manually
			//	collect them
				std::vector<Face*> vFaces;
				CollectFaces(vFaces, srcGrid, e, false);

			//	erase them
				for(std::vector<Face*>::iterator iter = vFaces.begin();
					iter != vFaces.end(); ++iter)
				{
					srcGrid.erase(*iter);
				}
			}

			if((!srcGrid.option_is_enabled(VOLOPT_AUTOGENERATE_EDGES)) &&
				(!srcGrid.option_is_enabled(VOLOPT_AUTOGENERATE_FACES)))
			{
			//	we have to erase them manually
			//	collect them
				std::vector<Volume*> vVolumes;
				CollectVolumes(vVolumes, srcGrid, e, false);

			//	erase them
				for(std::vector<Volume*>::iterator iter = vVolumes.begin();
					iter != vVolumes.end(); ++iter)
				{
					srcGrid.erase(*iter);
				}
			}
		//	erase the edge
			srcGrid.erase(e);
		}

	//	return the new vertex
		return newVertex;
	}

//	something went wrong in CreateEdgeSplitGeometry.
//	erase the new vertex and return NULL
	destGrid.erase(newVertex);
	return NULL;
}


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
		EdgeBase* e = *iter;
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

template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(EdgeBase* e, TVertexPositionAttachmentAccessor& aaPosVRT)
{
	typename TVertexPositionAttachmentAccessor::ValueType v;
//	init v with 0.
	VecSet(v, 0);

//	sum up
	VecAdd(v, aaPosVRT[e->vertex(0)], aaPosVRT[e->vertex(1)]);

//	average
	VecScale(v, v, 0.5);

	return v;
}

template <class TEdgeIterator>
void FixEdgeOrientation(Grid& grid, TEdgeIterator edgesBegin,
						TEdgeIterator edgesEnd)
{
	using namespace std;
	grid.begin_marking();
	
//	we'll mark all edges.
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
		grid.mark(*iter);
		
	stack<EdgeBase*> candidates;
	
//	iterate through the edges
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
	//	only marked edges have to be considered
		if(grid.is_marked(*iter)){
			grid.unmark(*iter);
			candidates.push(*iter);
			
			while(!candidates.empty()){
				EdgeBase* e = candidates.top();
				candidates.pop();
				
			//	iterate through all associated edges
				for(size_t i = 0; i < 2; ++i){
					for(Grid::AssociatedEdgeIterator assIter =
						grid.associated_edges_begin(e->vertex(i));
						assIter != grid.associated_edges_end(e->vertex(i)); ++assIter)
					{
						EdgeBase* ae = *assIter;
					//	fix orientation of marked associated edges.
					//	those edges are new candidates afterwards
						if(grid.is_marked(ae)){
						//	the local index of the vertex over which ae is connected to e has
						//	to be different from the local index of the vertex in e.
							if(GetVertexIndex(ae, e->vertex(i)) == (int)i){
							//	we have to flip the orientation
								grid.flip_orientation(ae);
							}
							
						//	prepare the new candidate
							grid.unmark(ae);
							candidates.push(ae);
						}
					}
				}
			}
		}
	}
	
	grid.end_marking();
}

template <class TEdgeIterator, class TAAPosVRT>
EdgeBase* FindShortestEdge(TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
							TAAPosVRT& aaPos)
{
//	if edgesBegin equals edgesEnd, then the list is empty and we can
//	immediately return NULL
	if(edgesBegin == edgesEnd)
		return NULL;

//	the first edge is the first candidate for the shortest edge.
//	We compare squares to avoid computation of the square root.
	EdgeBase* shortestEdge = *edgesBegin;
	number shortestLen = EdgeLengthSq(shortestEdge, aaPos);
	++edgesBegin;

	for(; edgesBegin != edgesEnd; ++edgesBegin){
		EdgeBase* curEdge = *edgesBegin;
		number curLen = EdgeLengthSq(curEdge, aaPos);
		if(curLen < shortestLen){
			shortestEdge = curEdge;
			shortestLen = curLen;
		}
	}

	return shortestEdge;
}

////////////////////////////////////////////////////////////////////////////////
template <class TEdgeIterator>
void RemoveDoubleEdges(Grid& grid, TEdgeIterator edgesBegin, TEdgeIterator edgesEnd)
{
//	iterate over all edges and check whether both associated vertices are already
//	marked. If this is the case and a marked edge exists between those vertices,
//	then the current edge is a double.

	grid.begin_marking();

//	the first is the double, the second the original
	std::vector<std::pair<EdgeBase*, EdgeBase*> > doubles;
	std::vector<EdgeBase*> edges;

	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter){
		EdgeBase* e = *iter;
		if(!grid.is_marked(e)){
		//	check whether both vertices are marked
			VertexBase* v0 = e->vertex(0);
			VertexBase* v1 = e->vertex(1);

			bool isDouble = false;

			if(grid.is_marked(v0) && grid.is_marked(v1)){
			//	a necessary condition is met. However not yet sufficient.
			//	find marked edge between v0 and v1.
				CollectAssociated(edges, grid, v0);
				for(size_t i = 0; i < edges.size(); ++i){
					EdgeBase* te = edges[i];
					if((te->vertex(0) == v1 || te->vertex(1) == v1)
						&& grid.is_marked(te))
					{
					//	e is a double
						isDouble = true;
						doubles.push_back(std::make_pair(e, te));
						break;
					}
				}
			}

		//	finally mark e and its vertices (every processed edge is marked).
			if(!isDouble){
				grid.mark(e);
				grid.mark(v0);
				grid.mark(v1);
			}
		}
	}

	grid.end_marking();

//	now erase all doubles
	for(size_t i = 0; i < doubles.size(); ++i){
	//	this allows listeners to take actions
		grid.objects_will_be_merged(doubles[i].second, doubles[i].second,
									doubles[i].first);
		grid.erase(doubles[i].first);
	}
}

}//	end of namespace

#endif

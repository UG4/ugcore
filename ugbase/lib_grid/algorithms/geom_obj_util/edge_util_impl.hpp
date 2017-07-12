/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__EDGE_UTIL_IMPL__
#define __H__LIB_GRID__EDGE_UTIL_IMPL__

//#include "edge_util.h"
#include <stack>
#include <queue>
#include <map>
#include "lib_grid/grid/grid_util.h"
#include "../orientation_util.h"
#include "vertex_util.h"
namespace ug
{

template <class face_iter_t>
void GetInnerEdgesOfFaceSoup(
			std::vector<Edge*>& edgesOut,
			Grid& g,
			face_iter_t facesBegin,
			face_iter_t facesEnd)
{
	using namespace std;
	map<Edge*, int>	m;

	Grid::edge_traits::secure_container edges;
	for(face_iter_t fiter = facesBegin; fiter != facesEnd; ++fiter)
	{
		g.associated_elements(edges, *fiter);
		for(size_t i = 0; i < edges.size(); ++i)
			++m[edges[i]];
	}

	for(map<Edge*, int>::iterator iter = m.begin(); iter != m.end(); ++iter)
	{
		if(iter->second > 1)
			edgesOut.push_back(iter->first);
	}
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
inline number EdgeLengthSq(const EdgeVertices* e, TAAPosVRT& aaPos)
{
	return VecDistanceSq(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
}

////////////////////////////////////////////////////////////////////////
template <class TAAPosVRT>
inline number EdgeLength(const EdgeVertices* e, TAAPosVRT& aaPos)
{
	return VecDistance(aaPos[e->vertex(0)], aaPos[e->vertex(1)]);
}

////////////////////////////////////////////////////////////////////////
//	SplitEdge
//	see edge_operations.h for detailed description
template<class TVertex>
TVertex* SplitEdge(Grid& grid, Edge* e, bool bConservative)
{
	return SplitEdge<TVertex>(grid, grid, e, NULL, bConservative);
}

////////////////////////////////////////////////////////////////////////
//	SplitEdge
//	see edge_operations.h for detailed description
template<class TVertex>
TVertex* SplitEdge(Grid& destGrid, Grid& srcGrid, Edge* e,
						AVertex* paAssociatedVertices,
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

template<class TVertexPositionAttachmentAccessor>
typename TVertexPositionAttachmentAccessor::ValueType
CalculateCenter(const Edge* e, TVertexPositionAttachmentAccessor& aaPosVRT)
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

////////////////////////////////////////////////////////////////////////
template<class TAAPosVRT, class TAAWeightVRT>
UG_API
typename TAAPosVRT::ValueType
CalculateCenter(const EdgeVertices* e, TAAPosVRT& aaPos, TAAWeightVRT& aaWeight)
{
	typename TAAPosVRT::ValueType v;
	typedef typename TAAWeightVRT::ValueType weight_t;

//	init v with 0.
	VecSet(v, 0);

//	sum up
	weight_t w0 = aaWeight[e->vertex(0)];
	weight_t w1 = aaWeight[e->vertex(1)];

	VecScaleAdd(v, w0, aaPos[e->vertex(0)], w1, aaPos[e->vertex(1)]);

	if((w0 + w1) != 0)
		VecScale(v, v, 1. / (number)(w0 + w1));

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
		
	stack<Edge*> candidates;
	
//	iterate through the edges
	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter)
	{
	//	only marked edges have to be considered
		if(grid.is_marked(*iter)){
			grid.unmark(*iter);
			candidates.push(*iter);
			
			while(!candidates.empty()){
				Edge* e = candidates.top();
				candidates.pop();
				
			//	iterate through all associated edges
				for(size_t i = 0; i < 2; ++i){
					for(Grid::AssociatedEdgeIterator assIter =
						grid.associated_edges_begin(e->vertex(i));
						assIter != grid.associated_edges_end(e->vertex(i)); ++assIter)
					{
						Edge* ae = *assIter;
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


template <class TEdgeIterator>
UG_API
void AdjustEdgeOrientationToFaceOrientation(Grid& grid, TEdgeIterator edgesBegin,
						   	   	   	   	    TEdgeIterator edgesEnd)
{
	using namespace std;

	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter){
		Face* nbr;
		int numNbrs = GetAssociatedFaces(&nbr, grid, *iter, 1);
		if(numNbrs != 1)
			continue;

		if(!EdgeOrientationMatches(*iter, nbr)){
			grid.flip_orientation(*iter);
		}
	}
}

template <class TEdgeIterator>
UG_API
void AdjustEdgeOrientationToFaceOrientation(Grid& grid, TEdgeIterator edgesBegin,
						   	   	   	   	    TEdgeIterator edgesEnd,
						   	   	   	   	    Grid::face_traits::callback considerFace)
{
	using namespace std;

	Grid::face_traits::secure_container faces;

	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter){
		Edge* e = *iter;
		grid.associated_elements(faces, e);

		size_t numConsidered = 0;
		Face* nbr = NULL;
		for(size_t i = 0; i < faces.size(); ++i){
			if(considerFace(faces[i])){
				++numConsidered;
				nbr = faces[i];
			}
		}

		if((numConsidered == 1) && (!EdgeOrientationMatches(e, nbr))){
			grid.flip_orientation(e);
		}
	}
}


template <class TEdgeIterator, class TAAPosVRT>
Edge* FindShortestEdge(TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
							TAAPosVRT& aaPos)
{
//	if edgesBegin equals edgesEnd, then the list is empty and we can
//	immediately return NULL
	if(edgesBegin == edgesEnd)
		return NULL;

//	the first edge is the first candidate for the shortest edge.
//	We compare squares to avoid computation of the square root.
	Edge* shortestEdge = *edgesBegin;
	number shortestLen = EdgeLengthSq(shortestEdge, aaPos);
	++edgesBegin;

	for(; edgesBegin != edgesEnd; ++edgesBegin){
		Edge* curEdge = *edgesBegin;
		number curLen = EdgeLengthSq(curEdge, aaPos);
		if(curLen < shortestLen){
			shortestEdge = curEdge;
			shortestLen = curLen;
		}
	}

	return shortestEdge;
}


template <class TEdgeIterator, class TAAPosVRT>
Edge* FindLongestEdge(TEdgeIterator edgesBegin, TEdgeIterator edgesEnd,
							TAAPosVRT& aaPos)
{
//	if edgesBegin equals edgesEnd, then the list is empty and we can
//	immediately return NULL
	if(edgesBegin == edgesEnd)
		return NULL;

//	the first edge is the first candidate for the shortest edge.
//	We compare squares to avoid computation of the square root.
	Edge* longestEdge = *edgesBegin;
	number longestLen = EdgeLengthSq(longestEdge, aaPos);
	++edgesBegin;

	for(; edgesBegin != edgesEnd; ++edgesBegin){
		Edge* curEdge = *edgesBegin;
		number curLen = EdgeLengthSq(curEdge, aaPos);
		if(curLen > longestLen){
			longestEdge = curEdge;
			longestLen = curLen;
		}
	}

	return longestEdge;
}
// ////////////////////////////////////////////////////////////////////////////////
// template <class TEdgeIterator>
// void RemoveDoubleEdges(Grid& grid, TEdgeIterator edgesBegin, TEdgeIterator edgesEnd)
// {
// //	iterate over all edges and check whether both associated vertices are already
// //	marked. If this is the case and a marked edge exists between those vertices,
// //	then the current edge is a double.

// 	grid.begin_marking();

// //	the first is the double, the second the original
// 	std::vector<std::pair<Edge*, Edge*> > doubles;
// 	std::vector<Edge*> edges;

// 	for(TEdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter){
// 		Edge* e = *iter;
// 		if(!grid.is_marked(e)){
// 		//	check whether both vertices are marked
// 			Vertex* v0 = e->vertex(0);
// 			Vertex* v1 = e->vertex(1);

// 			bool isDouble = false;

// 			if(grid.is_marked(v0) && grid.is_marked(v1)){
// 			//	a necessary condition is met. However not yet sufficient.
// 			//	find marked edge between v0 and v1.
// 				CollectAssociated(edges, grid, v0);
// 				for(size_t i = 0; i < edges.size(); ++i){
// 					Edge* te = edges[i];
// 					if((te->vertex(0) == v1 || te->vertex(1) == v1)
// 						&& grid.is_marked(te))
// 					{
// 					//	e is a double
// 						isDouble = true;
// 						doubles.push_back(std::make_pair(e, te));
// 						break;
// 					}
// 				}
// 			}

// 		//	finally mark e and its vertices (every processed edge is marked).
// 			if(!isDouble){
// 				grid.mark(e);
// 				grid.mark(v0);
// 				grid.mark(v1);
// 			}
// 		}
// 	}

// 	grid.end_marking();

// //	now erase all doubles
// 	for(size_t i = 0; i < doubles.size(); ++i){
// 	//	this allows listeners to take actions
// 		grid.objects_will_be_merged(doubles[i].second, doubles[i].second,
// 									doubles[i].first);
// 		grid.erase(doubles[i].first);
// 	}
// }

////////////////////////////////////////////////////////////////////////////////
template <class EdgeIterator, class TAAPos>
void MinimizeEdgeLength_SwapsOnly(Grid& grid, EdgeIterator edgesBegin,
								  EdgeIterator edgesEnd, TAAPos& aaPos)
{
	using namespace std;

//	helper to collect neighbors
	Face* nbrFaces[2];
	vector<Edge*> edges;

//	flipCandidates
	queue<Edge*> candidates;

//	sadly we can't use marking. Thats why we attach a simple byte to the edges,
//	which will tell whether an edge is already a candidate.
	AByte aIsCandidate;
	grid.attach_to_edges_dv(aIsCandidate, 0, false);
	Grid::AttachmentAccessor<Edge, AByte> aaIsCandidate(grid, aIsCandidate);

//	set up candidate array
	for(EdgeIterator iter = edgesBegin; iter != edgesEnd; ++iter){
		aaIsCandidate[*iter] = 1;
		candidates.push(*iter);
	}


	while(!candidates.empty()){
		Edge* e = candidates.front();
		candidates.pop();
		aaIsCandidate[e] = 0;

	//	we only perform swaps on regular manifolds.
		if(GetAssociatedFaces(nbrFaces, grid, e, 2) == 2){
		//	make sure that both neighbors are triangles
			if(nbrFaces[0]->num_vertices() != 3 || nbrFaces[1]->num_vertices() != 3)
				continue;

		//	check whether a swap would make the edge shorter.
			Vertex* conVrt0 = GetConnectedVertex(e, nbrFaces[0]);
			Vertex* conVrt1 = GetConnectedVertex(e, nbrFaces[1]);
			if(VertexDistanceSq(conVrt0, conVrt1, aaPos) < EdgeLengthSq(e, aaPos))
			{
			//	it'll be shorter
			//	now make sure that associated triangles won't flip
			//todo: add support for 2d position attachments
				vector3 n0, n1;
				CalculateNormal(n0, nbrFaces[0], aaPos);
				CalculateNormal(n1, nbrFaces[1], aaPos);
				number oldDot = VecDot(n0, n1);

				FaceDescriptor ntri;
				ntri.set_num_vertices(3);
				ntri.set_vertex(0, e->vertex(0));
				ntri.set_vertex(1, conVrt1);
				ntri.set_vertex(2, conVrt0);
				CalculateNormal(n0, &ntri, aaPos);

				ntri.set_vertex(0, e->vertex(1));
				ntri.set_vertex(1, conVrt0);
				ntri.set_vertex(2, conVrt1);
				CalculateNormal(n1, &ntri, aaPos);

				number newDot = VecDot(n0, n1);

			//	if both have the same sign, we're fine!
				if(oldDot * newDot < 0){
					continue;//	not fine!
				}

			//	ok - everything is fine. Now swap the edge
				e = SwapEdge(grid,  e);

				UG_ASSERT(e, "SwapEdge did not produce a new edge.");

			//	all edges of associated triangles are candidates again (except e)
				GetAssociatedFaces(nbrFaces, grid, e, 2);
				for(size_t i = 0; i < 2; ++i){
					CollectAssociated(edges, grid, nbrFaces[i]);
					for(size_t j = 0; j < edges.size(); ++j){
						if(edges[j] != e && (!aaIsCandidate[edges[j]])){
							candidates.push(edges[j]);
							aaIsCandidate[edges[j]] = 1;
						}
					}
				}
			}
		}
	}

	grid.detach_from_edges(aIsCandidate);
}

template <class vector_t, class TAAPos>
UG_API bool
ContainsPoint(const EdgeVertices* e, const vector_t& p, TAAPos aaPos)
{
	number center = (aaPos[e->vertex(0)].x() + aaPos[e->vertex(1)].x()) / 2.;
	number rad = fabs(aaPos[e->vertex(1)].x() - aaPos[e->vertex(0)].x()) / 2.;

	if(fabs(p.x() - center) <= rad)
		return true;
	return false;
}

}//	end of namespace

#endif

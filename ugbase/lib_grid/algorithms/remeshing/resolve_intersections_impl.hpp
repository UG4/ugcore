/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__resolve_intersections_impl__
#define __H__UG__resolve_intersections_impl__

#include <map>
#include <limits>
#include "resolve_intersections.h"
#include "common/math/misc/shapes.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/grid_generation/triangle_fill_sweep_line.h"
#include "lib_grid/algorithms/orientation_util.h"
#include "lib_grid/algorithms/remove_duplicates_util.h"
#include "lib_grid/algorithms/space_partitioning/lg_ntree.h"
#include "common/space_partitioning/ntree_traverser.h"

namespace ug{

template <class TAAPosVRT>
Vertex* ResolveVertexEdgeIntersection(Grid& grid, Vertex* v,
										   Edge* e, TAAPosVRT& aaPos,
										   number snapThreshold)
{
	typedef typename TAAPosVRT::ValueType vector_t;

	number snapThresholdSq = snapThreshold * snapThreshold;

//	make sure that the vertex is not an endpoint of e
	if(EdgeContains(e, v))
		return NULL;

//	we have to make sure that v and e are not connected by a face.
//	This could lead to infinite recursion
/*
	vector<Face*> faces;
	CollectFaces(faces, grid, e);
	for(size_t i = 0; i < faces.size(); ++i){
		if(FaceContains(faces[i], v))
			return NULL;
	}
*/
//	project the vertex to the line defined by the edge
	vector_t p;
	number t = DropAPerpendicular(p, aaPos[v], aaPos[e->vertex(0)],
								  aaPos[e->vertex(1)]);

	if((t >= 0) && (t <= 1.)){
		if(VecDistanceSq(p, aaPos[v]) < snapThresholdSq){
		//	to make sure that no double edges may occur, we'll use MergeVertices
			RegularVertex* nVrt = SplitEdge<RegularVertex>(grid, grid, e);
			aaPos[v] = p;
			MergeVertices(grid, v, nVrt);
			return v;
/*
		//	insert the vertex into the edge
			CreateEdgeSplitGeometry(grid, grid, e, v);
			grid.erase(e);
			return v;
*/
		}
	}
	return NULL;
}

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volumes, too.
 */
 template <class TAAPosVRT>
bool ResolveVertexFaceIntersection(Grid& grid, Vertex* v,
								   Face* f, TAAPosVRT& aaPos,
								   number snapThreshold,
								   std::vector<Face*>* pNewFacesOut)
{
	using namespace std;
	typedef typename TAAPosVRT::ValueType vector_t;

	if(pNewFacesOut)
		pNewFacesOut->clear();

	number snapThresholdSq = snapThreshold * snapThreshold;

//	make sure that the vertex is not a corner of f
	if(FaceContains(f, v))
		return false;

//	calculate the normal
	vector_t n;
	CalculateNormal(n, f, aaPos);

//	project the vertex to the plane defined by the face
	vector_t p;
	ProjectPointToPlane(p, aaPos[v], aaPos[f->vertex(0)], n);

//	check whether the distance is fine
	if(VecDistanceSq(p, aaPos[v]) < snapThresholdSq){
	//	now we have to check whether the projection lies in the face
		vector_t pi;
		bool refined = false;
		if(f->num_vertices() == 3){
			if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
										aaPos[f->vertex(2)], p, n))
			{
				refined = true;
			}
		}
		else if(f->num_vertices() == 4){
			if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
										aaPos[f->vertex(2)], p, n))
			{
				refined = true;
			}
			else if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(2)],
											aaPos[f->vertex(3)], p, n))
			{
				refined = true;
			}
		}

		if(!refined)
			return false;

	//	check whether the projection is too close to one of the corners.
	//	if this is the case, we'll set the position to that corner, so that
	//	a call to remove-doubles will resolve the intersection
		for(size_t i = 0; i < f->num_vertices(); ++i){
			if(VecDistanceSq(p, aaPos[f->vertex(i)]) < snapThresholdSq){
				// aaPos[v] = aaPos[f->vertex(i)];
				return false;
			}
		}

	//	adjust position
		//aaPos[v] = p;

	//	check whether the vertex is close to one of the edges of the face.
	//	if so, we'll only split the edge locally...
	//	this leads to T-junctions, which may be resolved later on
	//	through a call to ProjectVerticesToCloseEdges
		EdgeDescriptor ed;
		for(size_t i = 0; i < f->num_edges(); ++i){
			f->edge_desc(i, ed);
			vector_t pl;
			ProjectPointToLine(pl, p, aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]);
			if(VecDistanceSq(pl, p) < snapThresholdSq){
			//	we have to locally split the edge
				Vertex* newEdgeVrts[4] = {NULL, NULL, NULL, NULL};
				newEdgeVrts[i] = v;
				Vertex* newFaceVrt = NULL;
				vector<Face*> newFaces;
				if(f->refine(newFaces, &newFaceVrt, newEdgeVrts, NULL)){
					if(newFaceVrt)
						grid.register_element(newFaceVrt);
				//	split the edge to make sure that properties are passed on to children
					Edge* e = grid.get_edge(f, i);
					if(e){
						EdgeDescriptor ne1(ed.vertex(0), v);
						EdgeDescriptor ne2(v, ed.vertex(1));
						if(!grid.get_element(ne1))
							grid.create<RegularEdge>(ne1, e);
						if(!grid.get_element(ne2))
							grid.create<RegularEdge>(ne2, e);
					}
				//	insert the newly created faces into the grid
					for(size_t j = 0; j < newFaces.size(); ++j){
						if(!grid.get_face(*newFaces[j])){
							grid.register_element(newFaces[j], f);
							if(pNewFacesOut)
								pNewFacesOut->push_back(newFaces[j]);
						}
						else
							delete newFaces[j];
					}
				//	if f is the only connected face, we have to erase e
					if(e && NumAssociatedFaces(grid, e) == 1){
						grid.erase(f);
						grid.erase(e);
					}
					else
						grid.erase(f);
					return true;
				}
				else
					return false;

			}
		}

	//	if we reach this point we have to actually split the face.
	//	we do this by creating a new triangle connecting 'v' with each edge of f
		FaceDescriptor fd(3);
		for(size_t i = 0; i < f->num_edges(); ++i){
			f->edge_desc(i, ed);
			fd.set_vertex(0, ed.vertex(0));
			fd.set_vertex(1, ed.vertex(1));
			fd.set_vertex(2, v);
			if(!grid.get_element(fd)){
				Face* nf = *grid.create<Triangle>(TriangleDescriptor(fd.vertex(0),
												fd.vertex(1), fd.vertex(2)),
									  f);
				if(pNewFacesOut)
					pNewFacesOut->push_back(nf);
			}
		}
		grid.erase(f);
		return true;
	}

	return false;
}
// template <class TAAPosVRT>
// bool ResolveVertexFaceIntersection(Grid& grid, Vertex* v,
// 								   Face* f, TAAPosVRT& aaPos,
// 								   number snapThreshold,
// 								   std::vector<Face*>* pNewFacesOut)
// {
// 	using namespace std;
// 	typedef typename TAAPosVRT::ValueType vector_t;

// 	if(pNewFacesOut)
// 		pNewFacesOut->clear();

// 	number snapThresholdSq = snapThreshold * snapThreshold;

// //	make sure that the vertex is not a corner of f
// 	if(FaceContains(f, v))
// 		return false;

// //	calculate the normal
// 	vector_t n;
// 	CalculateNormal(n, f, aaPos);

// //	project the vertex to the plane defined by the face
// 	vector_t p;
// 	ProjectPointToPlane(p, aaPos[v], aaPos[f->vertex(0)], n);

// //	check whether the distance is fine
// 	if(VecDistanceSq(p, aaPos[v]) < snapThresholdSq){
// 		bool refined = false;
// 		vector<Face*> newFaces;
// 		Vertex* newFaceVrt = NULL;
// 		// Vertex* nVrt = NULL;
// 		vector_t pi;
// 	//	now we have to check whether the projection lies in the face
// 		if(f->num_vertices() == 3){
// 			if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
// 										aaPos[f->vertex(2)], p, n))
// 			{
// 				Vertex* newEdgeVrts[3] = {NULL, NULL, NULL};
// 				refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, v);
// 			}
// 		}
// 		else if(f->num_vertices() == 4){
// 			bool success = false;
// 			if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
// 										aaPos[f->vertex(2)], p, n))
// 			{
// 				success = true;
// 			}
// 			else if(RayTriangleIntersection(pi, aaPos[f->vertex(0)], aaPos[f->vertex(2)],
// 											aaPos[f->vertex(3)], p, n))
// 			{
// 				success = true;
// 			}

// 			if(success){
// 				Vertex* newEdgeVrts[4] = {NULL, NULL, NULL, NULL};
// 				refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, v);
// 			}
// 		}

// 		if(refined){
// 		//	adjust position
// 			aaPos[v] = pi;
// 		//	register the new faces and erase the old one
// 			for(size_t i = 0; i < newFaces.size(); ++i){
// 				if(!grid.get_face(*newFaces[i])){
// 					grid.register_element(newFaces[i], f);
// 					if(pNewFacesOut)
// 						pNewFacesOut->push_back(newFaces[i]);
// 				}
// 				else
// 					delete newFaces[i];
// 			}
// 			grid.erase(f);
// 			return true;
// 		}
// 	}

// 	return false;
// }

/**
 * This method does not resolve intersections between close, parallel edges or
 * between degenerate edges. You can treat such cases with
 * ReolveVertexEdgeIntersection.
 */
template <class TAAPosVRT>
Vertex* ResolveEdgeEdgeIntersection(Grid& grid, Edge* e1, Edge* e2,
										TAAPosVRT& aaPos, number snapThreshold)
{
	typedef typename TAAPosVRT::ValueType vector_t;

//	check whether one edge contains a vertex of another edge
	if(EdgeContains(e1, e2->vertex(0)) || EdgeContains(e1, e2->vertex(1)))
		return NULL;

	number snapThresholdSq = snapThreshold * snapThreshold;

	number t1, t2;
	if(LineLineProjection(t1, t2, aaPos[e1->vertex(0)], aaPos[e1->vertex(1)],
						  aaPos[e2->vertex(0)], aaPos[e2->vertex(1)]))
	{
	//	calculate the positions
		vector_t v1, v2;
		VecScaleAdd(v1, (1. - t1), aaPos[e1->vertex(0)], t1, aaPos[e1->vertex(1)]);
		VecScaleAdd(v2, (1. - t2), aaPos[e2->vertex(0)], t2, aaPos[e2->vertex(1)]);

	//	check whether the points are close to each other
		if(VecDistanceSq(v1, v2) < snapThresholdSq){
		//	calculate center
			vector_t p;
			VecScaleAdd(p, 0.5, v1, 0.5, v2);

		//	to make sure that no double edges may occur, we'll use MergeVertices
			RegularVertex* nVrt1 = SplitEdge<RegularVertex>(grid, grid, e1);
			RegularVertex* nVrt2 = SplitEdge<RegularVertex>(grid, grid, e2);
			aaPos[nVrt1] = p;
			MergeVertices(grid, nVrt1, nVrt2);

			return nVrt1;
		/*
		//	create a new vertex and split both edges using it
			RegularVertex* vrt = *grid.create<RegularVertex>();
			aaPos[vrt] = p;
			CreateEdgeSplitGeometry(grid, grid, e1, vrt);
			CreateEdgeSplitGeometry(grid, grid, e2, vrt);
			grid.erase(e1);
			grid.erase(e2);
		*/
		}
		else{
		/*
			LOG("distance check failed at: " << v1 << ", " << v2 << endl);
			UG_LOG("edges with vertices: " << aaPos[e1->vertex(0)] << aaPos[e1->vertex(1)] << endl);
			UG_LOG("                     " << aaPos[e2->vertex(0)] << aaPos[e2->vertex(1)]);
		*/
		}
	}
	return NULL;
}

/**
 * No support for volumes in the current version.
 * \todo Instead of manually refining the face, an external function SplitFace
 *		 should be used, which can take care of volume, too.
 */
template <class TAAPosVRT>
bool ResolveEdgeFaceIntersection(Grid& grid, Edge* e, Face* f,
								 TAAPosVRT& aaPos, number snapThreshold)
{
	using namespace std;
	typedef typename TAAPosVRT::ValueType vector_t;

//	check whether one edge contains a vertex of another edge
	if(FaceContains(f, e->vertex(0)) || FaceContains(f, e->vertex(1)))
		return false;

	vector_t dir;
	VecSubtract(dir, aaPos[e->vertex(1)], aaPos[e->vertex(0)]);

	vector_t p;
	number t1, t2, s;
	bool refined = false;
	vector<Face*> newFaces;
	Vertex* newFaceVrt = NULL;
	Vertex* vrt = NULL;
	if(f->num_vertices() == 3){
		if(RayTriangleIntersection(p, t1, t2, s, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
									aaPos[f->vertex(2)], aaPos[e->vertex(0)], dir))
		{
			if((s >= 0) && (s <= 1.)){
			//	split the face
				vrt = *grid.create<RegularVertex>();
				Vertex* newEdgeVrts[3] = {NULL, NULL, NULL};
				refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, vrt);
			}
		}
	}
	else if(f->num_vertices() == 4){
		bool intersecting = false;
		if(RayTriangleIntersection(p, t1, t2, s, aaPos[f->vertex(0)], aaPos[f->vertex(1)],
									aaPos[f->vertex(2)], aaPos[e->vertex(0)], dir))
		{
			intersecting = true;
		}
		else if(RayTriangleIntersection(p, t1, t2, s, aaPos[f->vertex(0)], aaPos[f->vertex(2)],
										aaPos[f->vertex(3)], aaPos[e->vertex(0)], dir))
		{
			intersecting = true;
		}

		if(intersecting && (s >= 0) && (s <= 1.))
		{
		//	split the face
			vrt = *grid.create<RegularVertex>();
			Vertex* newEdgeVrts[4] = {NULL, NULL, NULL, NULL};
			refined = f->refine(newFaces, &newFaceVrt, newEdgeVrts, vrt);
		}
	}

	if(refined && vrt){
	//	create a new vertex and adjust position
		aaPos[vrt] = p;

	//	register the new faces and erase the old one
		for(size_t i = 0; i < newFaces.size(); ++i)
			grid.register_element(newFaces[i], f);
		grid.erase(f);

	//	to make sure that no double edges may occur, we'll use MergeVertices
	//	and SplitEdge
		RegularVertex* nVrt = SplitEdge<RegularVertex>(grid, grid, e);
		MergeVertices(grid, vrt, nVrt);

/*
	//	split the edge with the new vertex and erase the old one
		CreateEdgeSplitGeometry(grid, grid, e, vrt);
		grid.erase(e);
*/

		return true;
	}

	return false;
}


///	sorts vertices along the specified ray
/**	inserts vertices into the specified multimap based on local coordinate of
 * their projection on the specified ray.
 * if clearContainer is specified as false (default is true), the vertices will
 * be inserted into the already existing sorted vertex set in vrtsOut.*/
template <class TVrtIter, class TAAPos>
void SpacialVertexSort(std::multimap<number, Vertex*>& vrtsOut,
						const TVrtIter vrtsBegin, const TVrtIter vrtsEnd,
						const typename TAAPos::ValueType& rayFrom,
						const typename TAAPos::ValueType& rayDir,
						TAAPos aaPos,
						bool clearContainer = true)
{
	typedef typename TAAPos::ValueType		vector_t;

	if(clearContainer)
		vrtsOut.clear();

	for(TVrtIter i_vrt = vrtsBegin; i_vrt != vrtsEnd; ++i_vrt){
		vector_t p;
		number t = ProjectPointToRay(p, aaPos[*i_vrt], rayFrom, rayDir);
		vrtsOut.insert(std::pair<number, Vertex*>(t, *i_vrt));
	}
}


/// Inserts the specified vertices on the given edge.
/**	The method assumes that the specified vertices do lie on the
 * specified edge. The specified vertices do not have to be sorted
 * in any specific order, since the method will sort them automatically.*/
template <class TVrtIter, class TAAPos>
void MultiEdgeSplit(Grid& grid, Edge* edge,
					TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					TAAPos aaPos)
{
	if(vrtsBegin == vrtsEnd)
		return;

	typedef typename TAAPos::ValueType		vector_t;
	vector_t dir;
	VecSubtract(dir, aaPos[edge->vertex(1)], aaPos[edge->vertex(0)]);

	std::multimap<number, Vertex*>	vrtMap;
	SpacialVertexSort(vrtMap, vrtsBegin, vrtsEnd, aaPos[edge->vertex(0)], dir, aaPos);

//	all associated faces have to be converted to triangles first
//	since we may alter faces, we'll collect them in a std::vector
	std::vector<Face*> faces;
	CollectFaces(faces, grid, edge);
	for(size_t i = 0; i < faces.size(); ++i){
		if(faces[i]->num_vertices() > 3){
			Quadrilateral* q = dynamic_cast<Quadrilateral*>(faces[i]);
			if(q)
				Triangulate(grid, q, &aaPos);
		}
	}

//	now collect associated faces again and find the connected vertices
//	also store for each associated face whether the orientation matches.
	CollectFaces(faces, grid, edge);
	std::vector<Vertex*>	conVrts(faces.size());
	std::vector<bool>		orientations(faces.size());
	for(size_t i = 0; i < faces.size(); ++i){
		orientations[i] = EdgeOrientationMatches(edge, faces[i]);
		conVrts[i] = GetConnectedVertex(edge, faces[i]);
		UG_COND_THROW(conVrts[i] == NULL,
					  " no connected vertex found to edge: " << ElementDebugInfo(grid, edge)
					  << " and face " << ElementDebugInfo(grid, faces[i]));
	}

	FaceDescriptor fd(3);
	TriangleDescriptor td;

	std::multimap<number, Vertex*>::iterator iter = vrtMap.begin();
	Vertex* curVrt = edge->vertex(0);
	Vertex* nextVrt = iter->second;
	while(nextVrt){
	//	create a new edge
		if(curVrt != nextVrt){
			EdgeDescriptor ed(curVrt, nextVrt);
			if(!grid.get_element(ed))
				grid.create<RegularEdge>(ed, edge);

		//	create new triangles for all connected vertices
			for(size_t i = 0; i < faces.size(); ++i){
				if(conVrts[i] == curVrt || conVrts[i] == nextVrt)
					continue;

				if(orientations[i]){
					fd.set_vertex(0, curVrt);
					fd.set_vertex(1, nextVrt);
					fd.set_vertex(2, conVrts[i]);
				}
				else{
					fd.set_vertex(0, nextVrt);
					fd.set_vertex(1, curVrt);
					fd.set_vertex(2, conVrts[i]);
				}
				if(!grid.get_element(fd)){
					grid.create<Triangle>(TriangleDescriptor(fd.vertex(0), fd.vertex(1), fd.vertex(2)),
										  faces[i]);
				}
			}
		}

	//	pick next vertex pair
		if(nextVrt == edge->vertex(1))
			break;
		else{
			curVrt = nextVrt;
			++iter;
			if(iter != vrtMap.end())
				nextVrt = iter->second;
			else
				nextVrt = edge->vertex(1);
		}
	}

//	erase the old edge and associated faces
	grid.erase(edge);
}


/**
 *	Projects vertices in elems onto close edges in elems.
 *	Though this method can be used to remove degenerated triangles,
 *	it is not guaranteed, that no degenerated triangles will remain
 *	(indeed, new degenerated triangles may be introduced).
 */
namespace impl{
namespace ProjectVerticesToCloseEdges{
	struct Record{
		Record() : closestEdge(-1), distanceSq(std::numeric_limits<number>::max()) {}
		int closestEdge;
		number distanceSq;
	};
}}

template <class TAPos>
bool ProjectVerticesToCloseEdges(Grid& grid,
								 GridObjectCollection elems,
								 TAPos& aPos,
								 number snapThreshold)
{
//	to speed things up we'll use an octree!
	typedef lg_ntree<3, 3, Vertex> octree_t;
	typedef octree_t::box_t	box_t;
	typedef typename TAPos::ValueType	vector_t;

	Grid::VertexAttachmentAccessor<TAPos> aaPos(grid, aPos);
	number snapThresholdSq = sq(snapThreshold);

	size_t numVrts = elems.num<Vertex>();
	size_t numEdges = elems.num<Edge>();

//	we'll only use an octree if there is a high number of vertices and edges...
	bool useOctree = (numVrts > 10 && numEdges > 10);
	octree_t octree(grid, aPos);
	if(useOctree)
		octree.create_tree(elems.begin<Vertex>(), elems.end<Vertex>());

//	we'll now iterate over all edges and search for close vertices.
	vector_t offset;
	VecSet(offset, snapThreshold);
	std::vector<Vertex*>	closeVrts;
	std::vector<Vertex*>	insertVrts;

	if(!useOctree){
		closeVrts.reserve(elems.num<Vertex>());
		for(VertexIterator iter = elems.begin<Vertex>();
			iter != elems.end<Vertex>(); ++iter)
		{
			closeVrts.push_back(*iter);
		}
	}

	using impl::ProjectVerticesToCloseEdges::Record;
	std::map<Vertex*, Record>	snapVrtMap;
	std::vector<Edge*> edges;
	edges.reserve(elems.num<Edge>());
	for(EdgeIterator eIter = elems.begin<Edge>(); eIter != elems.end<Edge>(); ++eIter)
		edges.push_back(*eIter);

//	find the closest edge for each vertex first and store it in snapVrtMap
	for(size_t iEdge = 0; iEdge < edges.size(); ++iEdge){
		Edge* e = edges[iEdge];

		vector_t corner[2];
		corner[0] = aaPos[e->vertex(0)];
		corner[1] = aaPos[e->vertex(1)];

	//	get all close vertices from the tree and do the projection
		if(useOctree){
			box_t bbox(corner, 2);
			bbox.min -= offset;
			bbox.max += offset;
			FindElementsInIntersectingNodes(closeVrts, octree, bbox);
		}

	//	find all vertices which have to be inserted by traversing closeVrts
		for(size_t i = 0; i < closeVrts.size(); ++i){
			Vertex* vrt = closeVrts[i];
			if(EdgeContains(e, vrt))
				continue;

			vector_t p;
			number s = ProjectPointToLine(p, aaPos[vrt], corner[0], corner[1]);
		//	we do an exact comparision here, since all other matches should be
		//	handled by a remove-doubles anyways
			if(s > 0 && s < 1.){
				const number distSq = VecDistanceSq(p, aaPos[vrt]);
				if(distSq <= snapThresholdSq){
					Record& rec = snapVrtMap[vrt];
					if(distSq < rec.distanceSq){
						rec.closestEdge = (int)iEdge;
						rec.distanceSq = distSq;
					}
				}
			}
		}
	}

//	split edges by finding associated vertices in snapVrtMap
	for(size_t iEdge = 0; iEdge < edges.size(); ++iEdge){
		Edge* e = edges[iEdge];

		vector_t corner[2];
		corner[0] = aaPos[e->vertex(0)];
		corner[1] = aaPos[e->vertex(1)];

	//	get all close vertices from the tree and do the projection
		if(useOctree){
			box_t bbox(corner, 2);
			bbox.min -= offset;
			bbox.max += offset;
			FindElementsInIntersectingNodes(closeVrts, octree, bbox);
		}

	//	find all vertices which have to be inserted by traversing closeVrts
		insertVrts.clear();
		for(size_t i = 0; i < closeVrts.size(); ++i){
			Vertex* vrt = closeVrts[i];
			if(EdgeContains(e, vrt))
				continue;

			std::map<Vertex*, Record>::iterator imap = snapVrtMap.find(vrt);
			if(imap != snapVrtMap.end()){
				Record& rec = imap->second;
				if(rec.closestEdge == (int)iEdge)
					insertVrts.push_back(vrt);
			}
		}

		if(!insertVrts.empty()){
			MultiEdgeSplit(grid, e, insertVrts.begin(), insertVrts.end(), aaPos);
		}
	}

	return true;
}

/**
 *	Projects vertices in elems onto close faces in elems.
 * TObjectCollection has to fulfill the interface of a GridObjectCollection.
 *
 * \note:	It is recommended to use a Grid ore a Selector as TObjectCollection,
 *			since a GridObjectCollection is static and it would thus not be
 *			possible to resolve all intersections.*/
template <class TObjectCollection, class TAPos>
bool ProjectVerticesToCloseFaces(Grid& grid,
								 TObjectCollection& elems,
								 TAPos& aPos,
								 number snapThreshold)
{
//	to speed things up we'll use an octree!
	typedef lg_ntree<3, 3, Vertex> octree_t;
	typedef octree_t::box_t	box_t;
	typedef typename TAPos::ValueType	vector_t;

	Grid::VertexAttachmentAccessor<TAPos> aaPos(grid, aPos);
	number snapThresholdSq = sq(snapThreshold);

	size_t numVrts = elems.template num<Vertex>();
	size_t numFaces = elems.template num<Face>();

//	we'll only use an octree if there is a high number of vertices and edges...
	bool useOctree = (numVrts > 10 && numFaces > 10);
	octree_t octree(grid, aPos);
	if(useOctree)
		octree.create_tree(elems.template begin<Vertex>(), elems.template end<Vertex>());

	vector_t offset;
	VecSet(offset, snapThreshold);
	
	std::vector<Vertex*>	closeVrts;
	std::vector<Vertex*>	insertVrts;

	if(!useOctree){
		closeVrts.reserve(elems.template num<Vertex>());
		for(VertexIterator iter = elems.template begin<Vertex>();
			iter != elems.template end<Vertex>(); ++iter)
		{
			closeVrts.push_back(*iter);
		}
	}

	std::queue<Face*> qFaces;

	for(FaceIterator iter = elems.template begin<Face>();
		iter != elems.template end<Face>(); ++iter)
	{
		qFaces.push(*iter);
	}

	std::vector<Face*> newFaces;
	// int maxIters = 1;
	while(!qFaces.empty()){
		// if(maxIters <= 0)
		// 	break;

		Face* f = qFaces.front();
		qFaces.pop();
		if(useOctree){
			box_t bbox(aaPos[f->vertex(0)], aaPos[f->vertex(0)]);
			for(size_t i = 1; i < f->num_vertices(); ++i){
				bbox = box_t(bbox, aaPos[f->vertex(i)]);
			}
			bbox.min -= offset;
			bbox.max += offset;
			FindElementsInIntersectingNodes(closeVrts, octree, bbox);
		}

		for(size_t i = 0; i < closeVrts.size(); ++i){
			Vertex* vrt = closeVrts[i];
			if(!FaceContains(f, vrt)){
			//	if the vertex is too close to a corner, we'll skip the insertion
				bool skip = false;
				for(size_t j = 0; j < f->num_vertices(); ++j){
					Vertex* corner = f->vertex(j);
					if(VecDistanceSq(aaPos[vrt], aaPos[corner]) <= snapThresholdSq){
						skip = true;
						break;
					}
				}

				if(skip)
					continue;

				if(ResolveVertexFaceIntersection(grid, vrt, f, aaPos, snapThresholdSq,
											  	 &newFaces))
				{
					// UG_LOG("resolved intersection of vertex at: " << aaPos[vrt] << std::endl);
					// --maxIters;
					for(size_t j = 0; j < newFaces.size(); ++j){
						qFaces.push(newFaces[j]);
					}
					break;
				}
			}
		}
	}

	return true;
}

/**THIS METHOD USES Grid::mark.
 * Intersects all edges in elems which are closer to each other
 * than snapThreshold.
 *
 * TObjectCollection has to fulfill the interface of a GridObjectCollection.
 * \todo:	speed up through octree (use MultiEdgeSplit)
 * \note:	It is recommended to use a Grid ore a Selector as TObjectCollection,
 *			since a GridObjectCollection is static and it would thus not be
 *			possible to resolve all intersections.*/
template <class TObjectCollection, class TAAPosVRT>
bool IntersectCloseEdges(Grid& grid,
						 TObjectCollection& elems,
						 TAAPosVRT& aaPos,
						 number snapThreshold)
{
//	we'll first mark all elements in elems to make sure that
//	only edges which were initially marked are intersected.
	grid.begin_marking();
	for(EdgeIterator iter = elems.template begin<Edge>();
		iter != elems.template end<Edge>(); ++iter)
	{
		grid.mark(*iter);
	}

//	perform edge/edge and edge/face intersections
	Grid::edge_traits::secure_container	edges;
//	if faces exist, more intersections may be necessary...
	const size_t maxIntersections = sq(elems.template num<Edge>());
	size_t numIntersections = 0;

	for(EdgeIterator mainIter = elems.template begin<Edge>();
		mainIter != elems.template end<Edge>();)
	{
		if(numIntersections > maxIntersections){
			UG_LOG("Intersection threshold reached in IntersectCloseEdges. "
				   " Not all intersections may have been resolved.\n");
			UG_LOG("  num-initial-edges: " << sqrt(maxIntersections) << "\n");
			grid.end_marking();
			return false;
		}

		Edge* e = *mainIter;
		++mainIter;

	//	only consider marked edges
		if(!grid.is_marked(e))
			continue;

	//	check all other edges up to e.
		for(EdgeIterator iter = elems.template begin<Edge>(); *iter != e;)
		{
			Edge* e2 = *iter;
			++iter;
			if(!grid.is_marked(e2))
				continue;

		//	if an intersection occured, we have to move on to the next edge in the queue,
		//	since the old edge no longer exists.
			Vertex* nVrt = ResolveEdgeEdgeIntersection(grid, e, e2, aaPos, snapThreshold);
			if(nVrt){
				++numIntersections;
			//	all edges connected to nVrt have to be marked, since they are new candidates!
				grid.associated_elements(edges, nVrt);
				for(size_t i = 0; i < edges.size(); ++i){
					grid.mark(edges[i]);
				}
				break;
			}
		}
	}
	grid.end_marking();
	return true;
}


///	returns the index of the first vertex closer to p than snapThreshold.
/**	returns -1 if nothing was found.*/
template <class TAAPosVRT>
int FindCloseVertexInArray(std::vector<Vertex*>& array,
							const typename TAAPosVRT::ValueType& p,
							TAAPosVRT& aaPos, number snapThreshold)
{
	number snapThrSq = snapThreshold * snapThreshold;
//	iterate over the array and check whether a vertex close to vrt already exists.
	for(size_t i = 0; i < array.size(); ++i){
		if(VecDistanceSq(aaPos[array[i]], p) < snapThrSq){
		//	we got one. return the index
			return (int)i;
		}
	}
	return -1;
}

////////////////////////////////////////////////////////////////////////
///	Intersects Coplanar Triangles
/**	fills a vector with the intersections on tri 1. Each pair of points
 * corresponds to an intersection-edge.
 * Please note that the returned list may contain multiple copies of a point.
 * \note	This method performs the full test and does not perform any optimizations
 *			like an initial bounding-box check. For best performance this method
 *			should thus only be called if the two triangles do most likely intersect.*/
inline bool IntersectCoplanarTriangles(std::vector<vector2>& edgesOut,
						const vector2& p00, const vector2& p01, const vector2& p02,
						const vector2& p10, const vector2& p11, const vector2& p12)
{
	// UG_LOG("intersecting...\n");

	edgesOut.clear();

	vector2 t0_ARR[] = {p00, p01, p02};
	vector2 t1_ARR[] = {p10, p11, p12};
	vector2* t0 = t0_ARR;
	vector2* t1 = t1_ARR;

//	to avoid rounding issues we search for the largest node-value and scale SMALL by that value
	number maxLenSq = 0;
	for(size_t i = 0; i < 3; ++i){
		maxLenSq = std::max(maxLenSq, VecLengthSq(t0[i]));
		maxLenSq = std::max(maxLenSq, VecLengthSq(t1[i]));
	}
	const number sml	 = SMALL * sqrt(maxLenSq);
	const number smlSq = sq(sml);
	// UG_LOG("sml: " << sml << ", smlSq: " << smlSq << "\n");

	// UG_LOG(" 2d-coords-t1: " << p00 << ", " << p01 << ", " << p02 << std::endl);
	// UG_LOG(" 2d-coords-t2: " << p10 << ", " << p11 << ", " << p12 << std::endl);
	// for(size_t i = 0; i < 3; ++i){
	// 	for(size_t j = 0; j < 3; ++j){
	// 		UG_LOG("VertexDistance " << i << "-" << j << ": " << VecDistance(t0[i], t1[j]) << std::endl);
	// 	}
	// }

//	first we check whether the two triangles are separable
//	note that t0 and t1 are switched twice, so that both are checked
	for(int i_tri = 0; i_tri < 2; ++i_tri){
		for(int i0 = 0; i0 < 3; ++i0){
			int i1 = (i0+1)%3;
			int i2 = (i0+2)%3;
			vector2 d;
			VecSubtract(d, t0[i1], t0[i0]);
			number refLen = VecLength(d);
			if(refLen < SMALL)
				continue;
			//vector2 n(d.y() / refLen, -d.x() / refLen);
			vector2 n(d.y(), -d.x());
			VecScale(n, n, (number) 1 / refLen);
			vector2 v;
			VecSubtract(v, t0[i2], t0[i0]);
			number refDot = VecDot(v, n);
			// UG_LOG("refDot: " << refDot << std::endl);
			bool separable = true;
			for(int j = 0; j < 3; ++j){
				VecSubtract(v, t1[j], t0[i0]);
				// UG_LOG(" first-pass:\n");
				// UG_LOG("  checkedDot: " << VecDot(v, n) << std::endl);
				// UG_LOG("  product: " << VecDot(v, n) * refDot << std::endl);
				// UG_LOG("  refLen*SMALL: " << refLen * SMALL << std::endl);
				//if(VecDot(v, n) * refDot > refLen * SMALL){
				if(VecDot(v, n) * refDot > sml){
					// UG_LOG("not separable\n");
				//	not separable along this edge
					separable = false;
					break;
				}
			}

			if(separable)
				return true;
		}
		std::swap(t0, t1);
	}

	bool isInside[3];
	int matchingCorner[] = {-1, -1, -1};//	matching corner in t0 for i-th point of t1

	for(int i = 0; i < 3; ++i){
		isInside[i] = PointIsInsideTriangle(t1[i], p00, p01, p02);
		//if(isInside[i]){
			for(int j = 0; j < 3; ++j){
				if(VecDistanceSq(t0[j], t1[i]) < smlSq){
					matchingCorner[i] = j;
					isInside[i] = true;
					// UG_LOG("matching corners: " << i << ", " << j << "\n");
					break;
				}
			}
		//}
	}

	for(int i0 = 0; i0 < 3; ++i0){
		int i1 = (i0+1)%3;
		if(isInside[i0]){
			if(isInside[i1]){
				if((matchingCorner[i0] == -1) || (matchingCorner[i1] == -1)){
				//	both points do lie in the triangle and thus the whole edge does
					edgesOut.push_back(t1[i0]);
					edgesOut.push_back(t1[i1]);
				}
			}
			else{
			//	an intersection with one of the other edges has to exist
				bool gotOne = false;
				for(int j0 = 0; j0 < 3; ++j0){
					int j1 = (j0+1)%3;
					vector2 vi;
					number s0, s1;
					if(LineLineIntersection2d(vi, s0, s1, t1[i0], t1[i1],
											  t0[j0], t0[j1], sml))
					{
					//	got our two points
						if(VecDistanceSq(t1[i0], vi) > smlSq){
							edgesOut.push_back(t1[i0]);
							edgesOut.push_back(vi);
							gotOne = true;
							break;
						}
					}
				}
				UG_COND_THROW(!gotOne && (matchingCorner[i0] == -1),
							  "An intersection has to exist if one point is "
							  "inside and one point is outside!\n"
							  << " checked line: " << t1[i0] << " - " << t1[i1]
						  	  << "\nchecked tri: " << t0[0] << ", " << t0[1] << ", " << t0[2]);
			}
		}
		else if(isInside[i1]){
			//	an intersection with one of the other edges has to exist
				bool gotOne = false;
				for(int j0 = 0; j0 < 3; ++j0){
					int j1 = (j0+1)%3;
					vector2 vi;
					number s0, s1;
					if(LineLineIntersection2d(vi, s0, s1, t1[i0], t1[i1],
											  t0[j0], t0[j1], sml))
					{
					//	got our two points
						if(VecDistanceSq(vi, t1[i1]) > smlSq){
							edgesOut.push_back(vi);
							edgesOut.push_back(t1[i1]);
							gotOne = true;
							break;
						}
					}
				}
				UG_COND_THROW(!gotOne && (matchingCorner[i1] == -1),
							  "An intersection has to exist if one point is "
							  "inside and one point is outside!\n"
							  << " checked line: " << t1[i0] << " - " << t1[i1]
						  	  << "\nchecked tri: " << t0[0] << ", " << t0[1] << ", " << t0[2]);
		}
		else{
		//	both points lie outside. An intersection may still exist. In this case
		//	however, the edge should intersect two other edges
			int numInts = 0;
			vector2 ints[2];
			for(int j0 = 0; j0 < 3; ++j0){
				int j1 = (j0+1)%3;
				number s0, s1;
				if(LineLineIntersection2d(ints[numInts], s0, s1, t1[i0], t1[i1],
										  t0[j0], t0[j1], sml))
				{
					++numInts;
					if(numInts > 1)
						break;
				}
			}
			UG_COND_THROW(numInts == 1, "Either 0 or 2 intersections have to exist!\n"
						  << " checked line: " << t1[i0] << " - " << t1[i1]
						  << "\nchecked tri: " << t0[0] << ", " << t0[1] << ", " << t0[2]);
			if(numInts == 2){
				edgesOut.push_back(ints[0]);
				edgesOut.push_back(ints[1]);
			}
		}
	}
	return !edgesOut.empty();
}

////////////////////////////////////////////////////////////////////////
///	Intersects Coplanar Triangles
/**	fills a vector with the intersections on tri 1. Each pair of points
 * corresponds to an intersection-edge.
 * Please note that the returned list may contain multiple copies of a point.
 * This method works even if the triangles are not perfectly coplanar. In this case
 * the projection of the triangle (p10,p11,p12) onto the triangle (p00,p01,p02) is
 * intersected.
 * Please note that the returned list may contain multiple copies of a point.
 * \note	This method performs the full test and does not perform any optimizations
 *			like an initial bounding-box check. For best performance this method
 *			should thus only be called if the two triangles do most likely intersect.*/
inline bool IntersectCoplanarTriangles(std::vector<vector3>& edgesOut,
						const vector3& p00, const vector3& p01, const vector3& p02,
						const vector3& p10, const vector3& p11, const vector3& p12)
{
	edgesOut.clear();

//	we have to project all points into the 2d space.
//	The transformation taken from the first triangle
	vector3 n;
	CalculateTriangleNormal(n, p00, p01, p02);
	matrix33 m;
	if(!ConstructOrthonormalSystem(m, n, 2)){
		UG_LOG("Construction of local orthonormal system failed...\n");
		return false;
	}
	
//	in order to make the transformation more robust, we'll rotate the triangles
//	around the center of tri-0
	vector3 c = p00;
	c += p01;
	c += p02;
	c *= (1./3.);
	
	// UG_LOG(" 3d-coords-t1: " << p00 << ", " << p01 << ", " << p02 << std::endl);
	// UG_LOG(" 3d-coords-t2: " << p10 << ", " << p11 << ", " << p12 << std::endl);
	// UG_LOG(" n: " << n << ", c: " << c << std::endl);

	vector3 v, vp0, vp1, vp2, vp3, vp4, vp5;
	VecSubtract(v, p00, c);
	TransposedMatVecMult(vp0, m, v);
	VecSubtract(v, p01, c);
	TransposedMatVecMult(vp1, m, v);
	VecSubtract(v, p02, c);
	TransposedMatVecMult(vp2, m, v);
	VecSubtract(v, p10, c);
	TransposedMatVecMult(vp3, m, v);
	VecSubtract(v, p11, c);
	TransposedMatVecMult(vp4, m, v);
	VecSubtract(v, p12, c);
	TransposedMatVecMult(vp5, m, v);

	std::vector<vector2> edges2d;
	if(IntersectCoplanarTriangles(edges2d,
				vector2(vp0.x(), vp0.y()), vector2(vp1.x(), vp1.y()), vector2(vp2.x(), vp2.y()),
				vector2(vp3.x(), vp3.y()), vector2(vp4.x(), vp4.y()), vector2(vp5.x(), vp5.y())))
	{
		edgesOut.resize(edges2d.size());
		for(size_t i = 0; i < edges2d.size(); ++i){
			vector3 v;
			VecCopy(v, edges2d[i], 0);
			MatVecMult(edgesOut[i], m, v);
			edgesOut[i] += c;
		}
		return true;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
///	Intersects Coplanar Triangles
/**	fills a vector with the intersections on tri 1. Each pair of points
 * corresponds to an intersection-edge.*/
template <class TAAPos>
bool IntersectCoplanarTriangles(std::vector<typename TAAPos::ValueType>& edgesOut,
								FaceVertices* tri1, FaceVertices* tri2,
								TAAPos aaPos)
{
	try{
		// UG_LOG("checking tris: " << CalculateCenter(tri1, aaPos)
		// 			<< ", " << CalculateCenter(tri2, aaPos) << std::endl);
		if((tri1->num_vertices() == 3) && (tri2->num_vertices() == 3)){

			return IntersectCoplanarTriangles(edgesOut,
					aaPos[tri1->vertex(0)], aaPos[tri1->vertex(1)], aaPos[tri1->vertex(2)],
					aaPos[tri2->vertex(0)], aaPos[tri2->vertex(1)], aaPos[tri2->vertex(2)]);
		}
	}
	catch(UGError& err){
		UG_THROW(err.get_msg() << "\n"
				 << "centers of involved triangles: " << CalculateCenter(tri1, aaPos)
		 		 << ", " << CalculateCenter(tri2, aaPos));
	}
	// UG_CATCH_THROW("centers of involved triangles: " << CalculateCenter(tri1, aaPos)
	// 				<< ", " << CalculateCenter(tri2, aaPos));
	return false;
}

template <class TAAPos>
Sphere<typename TAAPos::ValueType>
CalculateBoundingSphere(FaceVertices* face, TAAPos aaPos)
{
	Sphere<typename TAAPos::ValueType> s;
	s.center = CalculateCenter(face, aaPos);
	number maxDistSq = 0;
	for(size_t i = 0; i < face->num_vertices(); ++i){
		maxDistSq = std::max<number>(maxDistSq,
								VecDistanceSq(s.center, aaPos[face->vertex(i)]));
	}
	s.radius = sqrt(maxDistSq);
	return s;
}
////////////////////////////////////////////////////////////////////////
/**	This method uses Grid::mark
 */
template <class TAPos>
bool ResolveTriangleIntersections(Grid& grid, TriangleIterator trisBegin,
							  TriangleIterator trisEnd, number snapThreshold,
							  TAPos& aPos)
{
	using namespace std;
	// number snapThresholdSq = sq(snapThreshold);
//todo: add octree
//	we use a selector to select elements that shall be merged and
//	triangles that are to be processed and deleted.
	Selector sel(grid);
	sel.enable_autoselection(false);

	sel.select(trisBegin, trisEnd);
//	we first select all associated vertices and perform a merge on them
	SelectAssociatedVertices(sel, trisBegin, trisEnd);
	// SelectAssociatedEdges(sel, trisBegin, trisEnd);
	RemoveDoubles<3>(grid, sel.vertices_begin(), sel.vertices_end(),
					 aPos, snapThreshold);

////////////////////////////////
//	PERFORM AND RESOLVE TRIANGLE - TRIANGLE INTERSECTIONS
	Grid::VertexAttachmentAccessor<TAPos> aaPos(grid, aPos);

//	to speed things up we'll use an octree!
	typedef lg_ntree<3, 3, Face> octree_t;
	typedef octree_t::box_t	box_t;

	octree_t octree(grid, aPos);
	octree.create_tree(trisBegin, trisEnd);

//	clear edges and vertices from the selector. faces have to stay, since we will
//	operate on them now.
	sel.clear<Vertex>();
	sel.clear<Edge>();

//	enable selection inheritance, since we want new elements to be
//	selected in this selector
	sel.enable_selection_inheritance(true);

//	we need some attachments in order to store new vertices and edges for
//	each face.
	typedef Attachment<vector<Vertex*> >		AVrtVec;
	typedef Attachment<vector<pair<int, int> > >	AEdgeDescVec;
	AVrtVec aVrtVec;
	AEdgeDescVec aEdgeDescVec;
	grid.attach_to_faces(aVrtVec);
	grid.attach_to_faces(aEdgeDescVec);
	Grid::FaceAttachmentAccessor<AVrtVec> aaVrtVec(grid, aVrtVec);
	Grid::FaceAttachmentAccessor<AEdgeDescVec> aaEdgeDescVec(grid, aEdgeDescVec);

	std::vector<vector3> planarIntersections;

//	this vector will be used to find close triangles
	vector<Face*>	closeTris;
	Grid::vertex_traits::secure_container vrts;

//	iterate over all triangles and perform intersecion with other triangles
	for(TriangleIterator triIter1 = sel.begin<Triangle>();
		triIter1 != sel.end<Triangle>(); ++triIter1)
	{
		Face* t1 = *triIter1;

		box_t bbox;
		bbox.min = bbox.max = aaPos[t1->vertex(0)];
		for(size_t i = 1; i < t1->num_vertices(); ++i)
			bbox = box_t(bbox, aaPos[t1->vertex(i)]);
		bbox.min -= vector3(snapThreshold, snapThreshold, snapThreshold);
		bbox.max += vector3(snapThreshold, snapThreshold, snapThreshold);

	//	find close triangles
		FindElementsInIntersectingNodes(closeTris, octree, bbox);

	//	iterate over the rest of the triangles
		for(size_t i_close = 0; i_close < closeTris.size(); ++i_close){
			Face* t2 = closeTris[i_close];
			Face* t[2]; t[0] = t1; t[1] = t2;

		//	the attachments have a fixed order throughout the whole iteration...
		//	we want to make sure that each pair of triangles is checked only once.
		//	at the same time we avoid that t1 is intersected with t1...
			if(grid.get_attachment_data_index(t1) >= grid.get_attachment_data_index(t2))
				continue;

			Sphere<vector3> s1 = CalculateBoundingSphere(t1, aaPos);
			Sphere<vector3> s2 = CalculateBoundingSphere(t2, aaPos);
			if(VecDistance(s1.center, s2.center)
				> (s1.radius + s2.radius + snapThreshold + SMALL))
			{
				continue;
			}

		//todo:	Move both triangles closer to the origin to minimize rounding issues...
			
		//	perform normal comparision to handle coplanar triangles
			vector3 n1, n2;
			CalculateNormal(n1, t1, aaPos);
			CalculateNormal(n2, t2, aaPos);
			number d = VecDot(n1, n2);
			if(fabs(d) > 1. - SMALL){
			//	if the two triangles aren't in the same plane, there's nothing to do...
				if(DistancePointToPlane(aaPos[t2->vertex(0)], aaPos[t1->vertex(0)], n1) > snapThreshold)
					continue;
			//	perform coplanar triangle intersection for tri1 and tri2
			//	note that t1 and t2 are swapped twice - once at the end of each i_tri iteration.
				for(int i_tri = 0; i_tri < 2; ++i_tri){
					if(IntersectCoplanarTriangles(planarIntersections, t1, t2, aaPos)){
					//	we have to make sure that the corners of both triangles are
					//	contained in aVrtVec if an intersection occurs.
						for(int j_tri = 0; j_tri < 2; ++j_tri){
							vector<Vertex*>& vrts = aaVrtVec[t[j_tri]];
							if(vrts.empty()){
								for(size_t i = 0; i < t[j_tri]->num_vertices(); ++i)
									vrts.push_back(t[j_tri]->vertex(i));
							}
						}

						vector<Vertex*>& vrts = aaVrtVec[t1];

						for(size_t i = 0; i < planarIntersections.size(); i+=2){
							int inds[2];
							for(int j = 0; j < 2; ++j){
								inds[j] = FindCloseVertexInArray(vrts,
																 planarIntersections[i+j],
															   	 aaPos, snapThreshold);
								if(inds[j] == -1){
								//	check if the vertex is contained in the other tri...
									vector<Vertex*>& vrts2 = aaVrtVec[t2];
									int ind = FindCloseVertexInArray(vrts2,
														planarIntersections[i+j],
														aaPos, snapThreshold);
									if(ind != -1){
									//	insert the vertex into t1's list of vertices
										inds[j] = (int)vrts.size();
										vrts.push_back(vrts2[ind]);
									}
								}
								if(inds[j] == -1){
								//	we have to create a new vertex
									Vertex* vrt = *grid.create<RegularVertex>();
									aaPos[vrt] = planarIntersections[i+j];
									inds[j] = (int)vrts.size();
									vrts.push_back(vrt);
								}
							}
							if(inds[0] != inds[1])
								aaEdgeDescVec[t1].push_back(make_pair(inds[0], inds[1]));
						}
					}
				//	swap tris
					std::swap(t1, t2);
				}
				continue;
			}

		//	since the faces are not coplanar, we have to make sure
		//	that t1 and t2 do not share an edge (two vertices)
			size_t numShared = NumSharedVertices(grid, t1, t2);
			if(numShared > 1)
				continue;

			vector3 ip[2];
			if(TriangleTriangleIntersection(aaPos[t1->vertex(0)], aaPos[t1->vertex(1)],
											aaPos[t1->vertex(2)], aaPos[t2->vertex(0)],
											aaPos[t2->vertex(1)], aaPos[t2->vertex(2)],
											&ip[0], &ip[1], SMALL) == 1)
			{
				// UG_LOG("> DBG < Intersection points: " << ip[0] << ", " << ip[1] << endl);
			//	add an edge between the two points
			//	to avoid insertion of double points, we first check whether the point
			//	already exists in the triangle. Do this for both triangles.

			//	prepare both triangles.
			//todo: think about performance optimizations.
			//	insertion of corner points could be avoided by bloating the code a little.
			//	this could increase performance.
				for(size_t i_tri = 0; i_tri < 2; ++i_tri){
				//	If it is encountered for the first time,
				//	we'll add its corner-vertices to its list of vertices.
					Face* tri = t[i_tri];
					vector<Vertex*>& vrts = aaVrtVec[tri];
					if(vrts.empty()){
						for(size_t i = 0; i < tri->num_vertices(); ++i)
							vrts.push_back(tri->vertex(i));
					}

				// (a)	this version only checks for close vertices in the list
				//		of vertices of the currently processed triangle. This
				//		version will introduce more (and possibly very close)
				//		vertices, but should be more robust.
				// todo: only create vertices for each point once, for both tris!

					// int edgeInd[2] = {-1, -1};
					// for(size_t i_point = 0; i_point < 2; ++i_point){
					// 	const vector3& point = ip[i_point];
					// 	int closeInd = FindCloseVertexInArray(aaVrtVec[tri], point,
					// 								   		  aaPos, snapThreshold);
							
					// 	if(closeInd == -1){
					// 		Vertex* vrt = *grid.create<RegularVertex>();
					// 		aaPos[vrt] = point;
					// 		closeInd = (int)aaVrtVec[tri].size();
					// 		aaVrtVec[tri].push_back(vrt);
					// 	}

					// 	edgeInd[i_point] = closeInd;
					// }

					// if(edgeInd[0] != edgeInd[1])
					// 	aaEdgeDescVec[tri].push_back(
					// 			make_pair(edgeInd[0], edgeInd[1]));
				}

			//	(b) the version below also checks for existing vertices in other tris.
			//	however, this can lead to vertices which do not lie in the plane
			//	of an intersecting triangle. This in turn may lead to problems
			//	during triangulation.

				int inds1[2];
				int inds2[2];
				for(size_t i = 0; i < 2; ++i){
					int tind1 = FindCloseVertexInArray(aaVrtVec[t[0]], ip[i],
													   aaPos, snapThreshold);
					int tind2 = FindCloseVertexInArray(aaVrtVec[t[1]], ip[i],
													   aaPos, snapThreshold);
					if(tind1 == -1){
						if(tind2 == -1){
						//	we have to create a new vertex
							Vertex* vrt = *grid.create<RegularVertex>();
							aaPos[vrt] = ip[i];
							tind1 = (int)aaVrtVec[t[0]].size();
							tind2 = (int)aaVrtVec[t[1]].size();
							aaVrtVec[t[0]].push_back(vrt);
							aaVrtVec[t[1]].push_back(vrt);
						}
						else{
						//	the vertex already exists in t[1]
							tind1 = (int)aaVrtVec[t[0]].size();
							aaVrtVec[t[0]].push_back((aaVrtVec[t[1]])[tind2]);
						}
					}
					else if(tind2 == -1){
					//	the vertex already exists in t[0]
						tind2 = (int)aaVrtVec[t[1]].size();
						aaVrtVec[t[1]].push_back((aaVrtVec[t[0]])[tind1]);
					}

				//	ind1 now contains the index into the vertex array of t[0], at
				//	which a vertex with position ip[i] lies.
					inds1[i] = tind1;
					inds2[i] = tind2;
				}
			//	we found the indices of both endpoints and can now add an edge
			//	connecting both to the edgeDesc arrays of t[0] and t[1].
				if(inds1[0] != inds1[1])
					aaEdgeDescVec[t[0]].push_back(make_pair(inds1[0], inds1[1]));
				if(inds2[0] != inds2[1])
					aaEdgeDescVec[t[1]].push_back(make_pair(inds2[0], inds2[1]));
			}
		}
	}

//	all intersections have been resolved. Iterate over the triangles again and
//	create the new elements.
//	triangles that shall be deleted are pushed to vDelTris
	vector<Triangle*> vDelTris;
//	here we collect all vertices on which a merge has to be performed at the end
//	of the algorithm (vertices created through edge-edge intersections inside a triangle)
	vector<Vertex*> cutVertices;
	Grid tgrid(GRIDOPT_STANDARD_INTERCONNECTION);
	AInt aInt;
	AVertex aVrt;
	tgrid.attach_to_vertices(aPos);
	tgrid.attach_to_vertices(aInt);
	tgrid.attach_to_vertices_dv(aVrt, NULL);
	Grid::VertexAttachmentAccessor<TAPos> taaPos(tgrid, aPos);
	Grid::VertexAttachmentAccessor<AVertex> aaVrt(tgrid, aVrt);

//	holds vertices of tgrid, so that they are accessible by index.
	vector<Vertex*> tgridVrts;
//	we don't want the newly created faces to be selected
	sel.enable_selection_inheritance(false);

	std::vector<Vertex*> dblVrts;

	for(TriangleIterator triIter = sel.begin<Triangle>();
		triIter != sel.end<Triangle>(); ++triIter)
	{
		Triangle* t = *triIter;
	//	we only proceed if there are intersecion-edges at all
		if(aaVrtVec[t].size() > 3){
			tgrid.clear_geometry();
			tgridVrts.clear();

		//	copy vertices associated with t1 to tgrid
			vector<Vertex*>& vrts = aaVrtVec[t];
			for(size_t i = 0; i < vrts.size(); ++i){
				Vertex* vrt = *tgrid.create<RegularVertex>();
				aaVrt[vrt] = vrts[i];
				taaPos[vrt] = aaPos[vrts[i]];
				tgridVrts.push_back(vrt);
			}

		//	now create the edges. vertices are found by indexing tgridVrts
			vector<pair<int, int> >& edgeDescs = aaEdgeDescVec[t];

		//	tri edges
			tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[0], tgridVrts[1]));
			tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[1], tgridVrts[2]));
			tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[2], tgridVrts[0]));

		//	new edges
			for(size_t i = 0; i < edgeDescs.size(); ++i){
				tgrid.create<RegularEdge>(EdgeDescriptor(tgridVrts[edgeDescs[i].first],
												  tgridVrts[edgeDescs[i].second]));
			}

		//	due to construction there may be double-edges
			RemoveDuplicates(tgrid, tgrid.edges_begin(), tgrid.edges_end());

		//	we now have to resolve intersections between the edges
		//	first we'll try to snap vertices to edges
			ProjectVerticesToCloseEdges(tgrid, tgrid.get_grid_objects(),
										aPos, snapThreshold);

		//	remove doubles
		//	typically we only have a few vertices here...
		//todo:	use a special version of this algorithm which uses a kdtree if many
		//		vertices are present
			// dblVrts.clear();
			// for(VertexIterator iter = tgrid.begin<Vertex>();
			// 	iter != tgrid.end<Vertex>(); ++iter)
			// {
			// 	dblVrts.push_back(*iter);
			// }

			// const size_t numVrts = dblVrts.size();
			// for(size_t i = 0; i < numVrts; ++i){
			// 	Vertex* v0 = dblVrts[i];
			// 	if(!v0)	continue;
			// 	for(size_t j = i+1; j < numVrts; ++j){
			// 		Vertex* v1 = dblVrts[j];
			// 		if(!v1)	continue;
			// 		if(VecDistanceSq(aaPos[v0], aaPos[v1]) <= snapThresholdSq){
			// 			MergeVertices(tgrid, v0, v1);
			// 			dblVrts[j] = NULL;
			// 		}
			// 	}
			// }

			// //DEBUGGING BEGIN
			// {
			// 	UG_LOG("STORING TEST GRID\n");
			// 	static int dbgCounter = 0;
			// 	std::stringstream name;
			// 	name << "/home/sreiter/Desktop/_grids/failure-"
			// 		 << dbgCounter << ".ugx";
			// 	SaveGridToFile(tgrid, name.str().c_str(), aPos);

			// 	++dbgCounter;
			// }
			// //DEBUGGING END

		//	if there are only 3 vertices left, there's nothing to do...
			if(tgrid.num_vertices() <= 3)
				continue;
		//	now resolve edge/edge intersections
			//IntersectCloseEdges(tgrid, tgrid, taaPos, SMALL);
			IntersectCloseEdges(tgrid, tgrid, taaPos, snapThreshold);
			// if(!IntersectCloseEdges(tgrid, tgrid, taaPos, snapThreshold)){
			// 	static int dbgCounter = 0;
			// 	std::stringstream name;
			// 	name << "/home/sreiter/Desktop/failed_sweeplines/failed_intersect_close_edges_"
			// 		 << dbgCounter << ".ugx";
			// 	SaveGridToFile(tgrid, name.str().c_str(), aPos);

			// 	++dbgCounter;
			// }

		//	make sure that all vertices have an associated aaVrt
			for(VertexIterator viter = tgrid.vertices_begin();
				viter != tgrid.vertices_end(); ++viter)
			{
				if(!aaVrt[*viter]){
				//	since the vertex does not have an associated vertex in grid,
				//	it is clear that it has been created through an edge-edge cut.
				//	Associates of such vertices have to be merged later on.
					aaVrt[*viter] = *grid.create<RegularVertex>();
					aaPos[aaVrt[*viter]] = taaPos[*viter];
					cutVertices.push_back(aaVrt[*viter]);
				}
			}


		//	ok. Everything is prepared. We can now triangulate the grid.
			if(TriangleFill_SweepLine(tgrid, tgrid.edges_begin(), tgrid.edges_end(),
										aPos, aInt))
			{
			//	mark the triangle for deletion
				vDelTris.push_back(*triIter);

			//	add the triangles to the grid.
				for(TriangleIterator titer = tgrid.begin<Triangle>();
					titer != tgrid.end<Triangle>(); ++titer)
				{
					Triangle* ntri = *titer;
					FaceDescriptor fd(3);
					fd.set_vertex(0, aaVrt[ntri->vertex(0)]);
					fd.set_vertex(1, aaVrt[ntri->vertex(1)]);
					fd.set_vertex(2, aaVrt[ntri->vertex(2)]);

					if(!grid.get_element(fd)){
						grid.create<Triangle>(TriangleDescriptor(aaVrt[ntri->vertex(0)],
																aaVrt[ntri->vertex(1)],
																aaVrt[ntri->vertex(2)]),
											 *triIter);
					}
				}
			}
			else{

					UG_LOG("at:")
					size_t cnt = 0;
					for(VertexIterator iter = tgrid.begin<Vertex>();
						iter != tgrid.end<Vertex>(); ++iter, ++cnt)
					{
						if(cnt > 2)
							break;
						UG_LOG(" " << taaPos[*iter]);
					}
					UG_LOG(std::endl);

					// UG_THROW("ResolveTriangleIntersections failed!");

			// 	static int fileCounter = 1;
			// 	string filenamePrefix = "/home/sreiter/Desktop/failed_sweeplines/failed_sweepline_";
			// 	stringstream ss2d, ss3d;
			// 	ss2d << filenamePrefix << "2d_" << fileCounter << ".lgb";
			// 	ss3d << filenamePrefix << "3d_" << fileCounter << ".lgb";
			// 	++fileCounter;
			// 	UG_LOG("TriangleFill_SweepLine failed!\n");
			// 	SaveGridToFile(tgrid, ss3d.str().c_str(), aPos);
			// //	perform transformation to 2d and save that too.
			// 	std::vector<vector3> vrts;
			// 	for(VertexIterator iter = tgrid.vertices_begin();
			// 		iter != tgrid.vertices_end(); ++iter)
			// 	{
			// 		vrts.push_back(taaPos[*iter]);
			// 	}
			// 	std::vector<vector2> vrts2d(vrts.size());
			// 	TransformPointSetTo2D(&vrts2d.front(), &vrts.front(),
			// 						  vrts.size());

			// 	size_t counter = 0;
			// 	for(VertexIterator iter = tgrid.vertices_begin();
			// 		iter != tgrid.vertices_end(); ++iter, ++counter)
			// 	{
			// 		taaPos[*iter] = vector3(vrts2d[counter].x(), vrts2d[counter].y(), 0);
			// 	}

			// 	SaveGridToFile(tgrid, ss2d.str().c_str(), aPos);
			}
		}
	}

//	detach attachments (tgrid is deleted anyways)
	grid.detach_from_faces(aVrtVec);
	grid.detach_from_faces(aEdgeDescVec);

////////////////////////////////
//	GRID POSTPROCESS
//	before we merge vertices in cutVertices, we'll select all faces
//	in order to make sure that only valid faces will be deleted.
	sel.clear();
	sel.select(vDelTris.begin(), vDelTris.end());
	sel.select(cutVertices.begin(), cutVertices.end());

//	perform the merge (this has to be done on a selector.
//	  the current version of RemoveDoubles is a little restrictive
//	  in this regard.)
	if(!sel.empty<Vertex>()){
		RemoveDoubles<3>(grid, sel.vertices_begin(), sel.vertices_end(),
						 aPosition, snapThreshold);
	}

	sel.clear<Vertex>();
	sel.clear<Edge>();

//	finally delete all refined triangles and associated unused edges and vertices
	SelectInnerSelectionEdges(sel);
	SelectInnerSelectionVertices(sel);

	grid.erase(sel.begin<Face>(), sel.end<Face>());
	grid.erase(sel.begin<Edge>(), sel.end<Edge>());
	grid.erase(sel.begin<Vertex>(), sel.end<Vertex>());

	return true;
}

}// end of namespace

#endif

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

#ifndef __H__UG__manifold_smoothing__
#define __H__UG__manifold_smoothing__

#include "common/types.h"
#include "lib_grid/algorithms/geom_obj_util/face_util.h"
#include "../volume_calculation.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
template <typename TIterator, typename AAPosVRT>
void LaplacianSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations)
{
	Grid::edge_traits::secure_container edges;
	Grid::face_traits::secure_container faces;
	Grid::volume_traits::secure_container vols;

	bool gotFaces = grid.num<Face>() > 0;
	bool gotVols = grid.num<Volume>() > 0;

	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			Vertex* vrt = *iter;
			vector3 v;
			VecSet(v, 0);
			number weight = 0;

		//	calculate smoothing vector relative to neighbors
			grid.associated_elements(edges, vrt);
			for(size_t i = 0; i < edges.size(); ++i){
				// number w = EdgeLength(edges[i], aaPos);
				number w = 1;
				VecScaleAdd(v, 1, v, w, aaPos[GetConnectedVertex(edges[i], vrt)]);
				weight += w;				
			}

			if(gotFaces){
				grid.associated_elements(faces, vrt);
				for(size_t i = 0; i < faces.size(); ++i){
					Face* f = faces[i];
					// number a = CalculateVolume(f, aaPos);
					number a = 1;
					VecScaleAdd(v, 1, v, a,
						CalculateGridObjectCenter(
							grid.get_opposing_object(vrt, f),
							aaPos));

					weight += a;
				}
			}

			if(gotVols){
				grid.associated_elements(vols, vrt);
				for(size_t i = 0; i < vols.size(); ++i){
					Volume* vol = vols[i];
					// number a = CalculateVolume(vol, aaPos);
					number a = 1;
					VecScaleAdd(v, 1, v, a,
						CalculateGridObjectCenter(
							grid.get_opposing_object(vrt, vol),
							aaPos));

					weight += a;
				}
			}

			if(weight > 0){
				VecScale(v, v, 1. / weight);
				VecSubtract(v, v, aaPos[vrt]);
				VecScale(v, v, alpha);
				VecAdd(aaPos[vrt], aaPos[vrt], v);
			}
		}
	}
}

///	Smoothes vertices in a 2d-manifold in 3d-space be moving vertices in the
///	tangential plane, only.
/**	USES Grid::mark*/
template <typename TVrtIter, typename TAAPos3>
void TangentialSmoothSimple(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					  TAAPos3 aaPos, number alpha, size_t numIterations)
{
	Grid::traits<Face>::secure_container	faces;
	for(size_t iteration = 0; iteration < numIterations; ++iteration){
		for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter){
			Vertex* vrt = *iter;
			g.associated_elements(faces, vrt);
			if(faces.empty())
				continue;

		//	calculate the normal at the current vertex
			vector3 n(0, 0, 0);
			for(size_t i = 0; i < faces.size(); ++i){
				vector3 tn;
				CalculateNormal(tn, faces[i], aaPos);
				VecAdd(n, n, tn);
			}

		//	calculate the center of connected vertices
			vector3 c(0, 0, 0);
			size_t numConnected = 0;
			g.begin_marking();
			for(size_t i_face = 0; i_face < faces.size(); ++i_face){
				Face::ConstVertexArray vrts = faces[i_face]->vertices();
				const size_t numVrts = faces[i_face]->num_vertices();
				for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
					Vertex* v = vrts[i_vrt];
					if((v != vrt) && (!g.is_marked(v))){
						g.mark(v);
						VecAdd(c, c, aaPos[v]);
						++numConnected;
					}
				}
			}
			g.end_marking();
			if(numConnected == 0)
				continue;

			VecScale(c, c, 1. / (number)numConnected);

		//	project the center of connected vertices to the plane through vrt
			vector3 cp;
			ProjectPointToPlane(cp, c, aaPos[vrt], n);

		//	move the vertex
			VecScaleAdd(aaPos[vrt], (1. - alpha), aaPos[vrt], alpha, cp);
		}
	}
}

///	Smoothes vertices in a 2d-manifold in 3d-space be moving vertices in the
///	tangential plane, only.
/**	USES Grid::mark
 * Computes an offset vector for each vertex and smoothes this offset by averaging
 * with offsets of adjacent vertices. The offset is then projected back into the
 * tangential plane and each vertex is relocated using those smoothed offsets.*/
template <typename TVrtIter, typename TAAPos3>
void TangentialSmooth(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
					  TAAPos3 aaPos, number alpha, size_t numIterations)
{
//	attach an integer index to all vertices
	AInt aInt;
	g.attach_to_vertices_dv(aInt, -1);
	Grid::VertexAttachmentAccessor<AInt> aaInt(g, aInt);

	int counter = 0;
	std::vector<vector3*> vrtPos;
	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter){
		aaInt[*iter] = counter++;
		vrtPos.push_back(&aaPos[*iter]);
	}

	std::vector<int> adjVrtInds;
	std::vector<int> adjVrtOffsets;
	std::vector<Face*> adjFaces;
	std::vector<int> adjFaceOffsets;

	adjVrtInds.reserve(counter * 7);//just a guess
	adjVrtOffsets.reserve(counter + 1);
	vrtPos.reserve(counter);
	adjFaces.reserve(counter * 7);//just a guess
	adjFaceOffsets.reserve(counter + 1);

//	build adjacency graph
	Grid::traits<Face>::secure_container	faces;
	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter){
		Vertex* vrt = *iter;
		adjVrtOffsets.push_back(adjVrtInds.size());
		adjFaceOffsets.push_back(adjFaces.size());

		g.associated_elements(faces, vrt);
		if(faces.empty())
			continue;

	//	create adjacencies
		g.begin_marking();
		for(size_t i_face = 0; i_face < faces.size(); ++i_face){
			adjFaces.push_back(faces[i_face]);
			Face::ConstVertexArray vrts = faces[i_face]->vertices();
			const size_t numVrts = faces[i_face]->num_vertices();
			for(size_t i_vrt = 0; i_vrt < numVrts; ++i_vrt){
				Vertex* v = vrts[i_vrt];
				if((v != vrt) && (!g.is_marked(v))){
					g.mark(v);
					if(aaInt[v] == -1){
					//	positions of vertices which are not to be smoothed but which are
					//	adjacent to smoothed ones, are inserted on the fly.
						aaInt[v] = (int)vrtPos.size();
						vrtPos.push_back(&aaPos[v]);
					}
					adjVrtInds.push_back(aaInt[v]);
				}
			}
		}
		g.end_marking();
	}

	adjVrtOffsets.push_back(adjVrtInds.size());
	adjFaceOffsets.push_back(adjFaces.size());
	g.detach_from_vertices(aInt);

	std::vector<vector3> offsets(vrtPos.size());
	std::vector<vector3> newOffsets(vrtPos.size());

	for(size_t iteration = 0; iteration < numIterations; ++iteration){
	//	calculate the initial offset vectors for each vertex
		for(size_t i_vrt = 0; i_vrt + 1 < adjVrtOffsets.size(); ++i_vrt){
		//	calculate plane normal
			int adjFacesBegin = adjFaceOffsets[i_vrt];
			int adjFacesEnd = adjFaceOffsets[i_vrt+1];

			if(adjFacesBegin == adjFacesEnd){
				offsets[i_vrt] = vector3(0, 0, 0);
				continue;
			}

			vector3 n(0, 0, 0);
			for(int i_adj = adjFacesBegin; i_adj < adjFacesEnd; ++i_adj){
				vector3 tn;
				CalculateNormal(tn, adjFaces[i_adj], aaPos);
				VecAdd(n, n, tn);
			}

		//	calculate center of adjacent vertices
			int adjVrtsBegin = adjVrtOffsets[i_vrt];
			int adjVrtsEnd = adjVrtOffsets[i_vrt+1];

			if(adjVrtsBegin == adjVrtsEnd){
				offsets[i_vrt] = vector3(0, 0, 0);
				continue;
			}

			vector3 c(0, 0, 0);
			for(int i_adj = adjVrtsBegin; i_adj < adjVrtsEnd; ++i_adj){
				int adjVrt = adjVrtInds[i_adj];
				VecAdd(c, c, *vrtPos[adjVrt]);
			}
			VecScale(c, c, 1. / number(adjVrtsEnd - adjVrtsBegin));

		//	project center to plane and calculate the offset
			vector3 cp;
			ProjectPointToPlane(cp, c, *vrtPos[i_vrt], n);
			VecSubtract(offsets[i_vrt], cp, *vrtPos[i_vrt]);
		}

	//	now smooth the offset vectors by repeatedly taking the average of
	//	neighbored vectors
		const size_t numOffsetSmoothSteps = 1;
		for(size_t i_step = 0; i_step < numOffsetSmoothSteps; ++i_step){
			for(size_t i_vrt = 0; i_vrt + 1 < adjVrtOffsets.size(); ++i_vrt){
				int adjVrtsBegin = adjVrtOffsets[i_vrt];
				int adjVrtsEnd = adjVrtOffsets[i_vrt+1];

				if(adjVrtsBegin == adjVrtsEnd)
					continue;

				int numAdj = adjVrtsEnd - adjVrtsBegin;
				//vector3 o = offsets[i_vrt];
				vector3 o;
				VecScale(o, offsets[i_vrt], numAdj);
				for(int i_adj = adjVrtsBegin; i_adj < adjVrtsEnd; ++i_adj){
					int adjVrt = adjVrtInds[i_adj];
					VecAdd(o, o, offsets[adjVrt]);
				}
				VecScale(newOffsets[i_vrt], o, 1. / number(2 * numAdj));

			//	restrict offset to plane
				int adjFacesBegin = adjFaceOffsets[i_vrt];
				int adjFacesEnd = adjFaceOffsets[i_vrt+1];

				if(adjFacesBegin == adjFacesEnd){
					offsets[i_vrt] = vector3(0, 0, 0);
					continue;
				}

				vector3 n(0, 0, 0);
				for(int i_adj = adjFacesBegin; i_adj < adjFacesEnd; ++i_adj){
					vector3 tn;
					CalculateNormal(tn, adjFaces[i_adj], aaPos);
					VecAdd(n, n, tn);
				}

				VecNormalize(n, n);
				number dot = VecDot(n, newOffsets[i_vrt]);
				VecScale(n, n, dot);
				VecSubtract(newOffsets[i_vrt], newOffsets[i_vrt], n);
			}
			offsets.swap(newOffsets);
		}

	//	and now move the vertices by their smoothed offset vectors
		for(size_t i_vrt = 0; i_vrt + 1 < adjVrtOffsets.size(); ++i_vrt){
			VecScaleAdd(*vrtPos[i_vrt], 1, *vrtPos[i_vrt], alpha, offsets[i_vrt]);
		}
	}
}


////////////////////////////////////////////////////////////////////////
/** vertices which will not be smoothed get a special weight when being considered
 * during smoothing of neighbored vertices.
 *
 * Make sure that cbSmoothVertex returns true for exactly all those vertices,
 * which in the range of the specified iterators.
 */
template <typename TIterator, typename AAPosVRT>
void WeightedEdgeSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations,
					Grid::vertex_traits::callback cbSmoothVertex)
{
	using vector_t = typename AAPosVRT::ValueType;

	Grid::edge_traits::secure_container edges;

	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			Vertex* vrt = *iter;
			vector_t vrtPos = aaPos[vrt];
			vector_t avDir;
			VecSet(avDir, 0);
			number weight = 0;

		//	calculate smoothing vector relative to neighbors
			grid.associated_elements(edges, vrt);
			number numNonSmooth = 0;
			for(size_t i = 0; i < edges.size(); ++i){
				Vertex* connVrt = GetConnectedVertex(edges[i], vrt);
				if(!cbSmoothVertex(connVrt))
					numNonSmooth += 1;
			}

			number nonSmoothWeight = 1. / std::max<number>(1, numNonSmooth);

			for(size_t i = 0; i < edges.size(); ++i){
				Vertex* connVrt = GetConnectedVertex(edges[i], vrt);
				number w = 1;
				if(!cbSmoothVertex(connVrt))
					w = nonSmoothWeight;

				vector_t dir;
				VecSubtract(dir, aaPos[connVrt], vrtPos);
				w *= VecLengthSq(dir);
				dir *= w;
				VecAdd(avDir, avDir, dir);
				weight += w;
			}

			if(weight > 0){
				avDir *= alpha / weight;
				VecAdd(aaPos[vrt], vrtPos, avDir);
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////
/** vertices which will not be smoothed get a special weight when being considered
 * during smoothing of neighbored vertices.
 *
 * Make sure that cbSmoothVertex returns true for exactly all those vertices,
 * which in the range of the specified iterators.
 */
template <typename TIterator, typename AAPosVRT>
void WeightedFaceSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, int numIterations,
					Grid::vertex_traits::callback cbSmoothVertex)
{
	using vector_t = typename AAPosVRT::ValueType;

	Grid::face_traits::secure_container faces;

	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			Vertex* vrt = *iter;
			vector_t vrtPos = aaPos[vrt];
			vector_t avDir;
			VecSet(avDir, 0);
			number weight = 0;

		//	calculate smoothing vector relative to neighbors
			grid.associated_elements(faces, vrt);

			for(size_t i = 0; i < faces.size(); ++i){
				Face* f = faces[i];
				if(f->num_vertices() != 3)
					continue;

				number w = 1;
				Edge* e = GetConnectedEdge(grid, vrt, f);
				Vertex* v0 = e->vertex(0);
				Vertex* v1 = e->vertex(1);
				vector_t edgeCenter;
				if(cbSmoothVertex(v0)){
					if(cbSmoothVertex(v1))
						VecScaleAdd(edgeCenter, 0.5, aaPos[v0], 0.5, aaPos[v1]);
					else
						VecScaleAdd(edgeCenter, 0.75, aaPos[v0], 0.25, aaPos[v1]);
				}
				else if(cbSmoothVertex(v1))
					VecScaleAdd(edgeCenter, 0.25, aaPos[v0], 0.75, aaPos[v1]);
				else{
					VecScaleAdd(edgeCenter, 0.5, aaPos[v0], 0.5, aaPos[v1]);
					w = 0.5;
				}

				vector_t dir;
				VecSubtract(dir, edgeCenter, vrtPos);
				w *= CalculateVolume(f, aaPos);
				dir *= w;
				VecAdd(avDir, avDir, dir);
				weight += w;
			}

			if(weight > 0){
				avDir *= alpha / weight;
				VecAdd(aaPos[vrt], vrtPos, avDir);
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////
/** The higher the dot-product between an outgoing edge and the vertex normal,
 * the higher the influence of that edge during smoothing of that vertex.
 */
template <typename TIterator, typename AAPosVRT>
void WeightedNormalSmooth(Grid& grid, TIterator vrtsBegin,
						  TIterator vrtsEnd, AAPosVRT& aaPos,
						  number alpha, number dotThreshold,
						  int numIterations)
{
	using std::min;
	using std::max;

	using vector_t = typename AAPosVRT::ValueType ;

	Grid::edge_traits::secure_container edges;

	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			Vertex* vrt = *iter;
			vector_t vrtPos = aaPos[vrt];
			vector_t vrtNorm;
			CalculateVertexNormal(vrtNorm, grid, vrt, aaPos);
			
			vector_t avDir;
			VecSet(avDir, 0);

			grid.associated_elements(edges, vrt);
			number nbrs = 0;
			for(size_t i = 0; i < edges.size(); ++i){
				Vertex* connVrt = GetConnectedVertex(edges[i], vrt);
				vector_t dir = aaPos[connVrt];
				dir -= vrtPos;
				vector_t ndir;
				VecNormalize(ndir, dir);
				number dot = VecDot(vrtNorm, ndir);
				if(dot >= dotThreshold){
					dir *= 0.5 * (dot + 1);
					avDir += dir;
					++nbrs;
				}
			}

			if(nbrs > 0){
				avDir *= alpha / nbrs;
				aaPos[vrt] += avDir;
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////
/** Linearizes the geometry along the gradient (steepest descent) with respect to the given up-vector*/
template <typename TIterator, typename AAPosVRT>
void SlopeSmooth(Grid& grid, TIterator vrtsBegin,
					TIterator vrtsEnd, AAPosVRT& aaPos,
					number alpha, const vector3& up,
					int numIterations)
{
	using std::min;
	using std::max;

	using vector_t = typename AAPosVRT::ValueType ;

	vector3 upN;
	VecNormalize(upN, up);

//	we'll approximate the gradient by finding the path between the intersections
//	of connected edges and the plane through the center-vertex, spanned by the normal-
//	and up-vectors
	Grid::edge_traits::secure_container edges;
	Grid::face_traits::secure_container faces;

	for(int iteration = 0; iteration < numIterations; ++iteration){
	//	iterate through all vertices
		for(TIterator iter = vrtsBegin; iter != vrtsEnd; ++iter){
		//	smooth each one
			Vertex* vrt = *iter;
			vector_t vrtPos = aaPos[vrt];
			vector_t vrtNorm;

		//	check if the vertex is a boundary vertex. If so, we'll only consider
		//	the connected boundary vertices
			if(LiesOnBoundary(grid, vrt)){
				grid.associated_elements(edges, vrt);

				vector3 target(0, 0, 0);
				number numBndEdges = 0;

				for(size_t i = 0; i < edges.size(); ++i){
					if(IsBoundaryEdge2D(grid, edges[i])){
						target += aaPos[GetConnectedVertex(edges[i], vrt)];
						++numBndEdges;
					}
				}
				
				if(numBndEdges > 0)
					VecScaleAdd(aaPos[vrt], alpha/numBndEdges, target, 1. - alpha, vrtPos);
			}
			else{
				CalculateVertexNormal(vrtNorm, grid, vrt, aaPos);
				
			//todo: one could think about performing laplacian smoothing in this case
				if(fabs(VecDot(vrtNorm, upN)) > 1. - SMALL)
					continue;

				vector3 planeNormal;
				VecCross(planeNormal, vrtNorm, upN);

				vector3	inters[2];
				int numInters = 0;

				grid.associated_elements(faces, vrt);
				for(size_t i = 0; i < faces.size(); ++i){
					Edge* e = GetConnectedEdge(grid, vrt, faces[i]);
					vector3 v0 = aaPos[e->vertex(0)];
					vector3 dir = aaPos[e->vertex(1)];
					dir -= v0;

					vector3 vi;
					number t;
					if(RayPlaneIntersection(vi, t, v0, dir, vrtPos, planeNormal)){
						if(t >= -SMALL && t <= 1 + SMALL){
							if((numInters > 1) && 
							   (VecDistanceSq(vi, inters[0]) > VecDistanceSq(inters[1], inters[0])))
							{
								inters[1] = vi;
							}
							else{
								inters[numInters] = vi;
								++numInters;
							}
						}
					}
				}

				if(numInters == 2){
					vector3 target;
					VecScaleAdd(target, 0.5, inters[0], 0.5, inters[1]);
					VecScaleAdd(aaPos[vrt], alpha, target, 1. - alpha, vrtPos);
				}
			}
		}
	}
}


}// end of namespace

#endif

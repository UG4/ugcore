/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_simplification
#define __H__UG_simplification

#include "common/math/ugmath.h"
#include "lib_grid/algorithms/debug_util.h"
#include "lib_grid/algorithms/polychain_util.h"
#include "lib_grid/grid/grid.h"

namespace ug{

inline bool Valence3VertexIsObsolete(const vector3& p, const vector3& c0,
							  const vector3& c1, const vector3& c2,
							  number maxSquaredHeightToBaseAreaRatio)
{
	vector3 proj, n;
	CalculateTriangleNormalNoNormalize(n, c0, c1, c2);
	if(RayTriangleIntersection(proj, c0, c1, c2, p, n)){
		number hsq = VecDistanceSq(proj, p);
		number area = TriangleArea(c0, c1, c2);

		if(area > 0)
			return (hsq / area) < maxSquaredHeightToBaseAreaRatio;
	}
	return false;
}

/** diagonal is always assumed to be between c0 and c2*/
inline bool Valence4VertexIsObsolete(const vector3& p, const vector3& c0,
									  const vector3& c1, const vector3& c2,
									  const vector3& c3,
									  number maxSquaredHeightToBaseAreaRatio)
{
	vector3 proj, n, n1, n2;
	CalculateTriangleNormalNoNormalize(n1, c0, c1, c2);
	CalculateTriangleNormalNoNormalize(n2, c0, c2, c3);

	VecAdd(n, n1, n2);

	if(RayTriangleIntersection(proj, c0, c1, c2, p, n)){
		number hsq = VecDistanceSq(proj, p);
		number area = TriangleArea(c0, c1, c2) + TriangleArea(c0, c2, c3);

		if(area > 0)
			return (hsq / area) < maxSquaredHeightToBaseAreaRatio;
	}

	if(RayTriangleIntersection(proj, c0, c2, c3, p, n)){
		number hsq = VecDistanceSq(proj, p);
		number area = TriangleArea(c0, c1, c2) + TriangleArea(c0, c2, c3);

		if(area > 0)
			return (hsq / area) < maxSquaredHeightToBaseAreaRatio;
	}

	return false;
}

///	returns true if the normals of the resulting triangles point into the same direction
/**	the quad is always assumed to be split along the diagonal c0, c2*/
inline bool QuadDiagonalIsValid(const vector3& c0, const vector3& c1,
								const vector3& c2, const vector3& c3)
{
	vector3 n0, n1;
	CalculateTriangleNormalNoNormalize(n0, c0, c1, c2);
	CalculateTriangleNormalNoNormalize(n1, c0, c2, c3);
	return VecDot(n0, n1) > SMALL;
}


template <class TVrtIter, class TAAPos>
void ReplaceValence3Vertices(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
							 number maxSquaredHeightToBaseAreaRatio, TAAPos aaPos)
{
	using namespace std;

	Grid::face_traits::secure_container	faces;
	Vertex* conVrts[3];

	std::vector<Vertex*> candidates;

	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd;)
	{
		Vertex* vrt = *iter;
		++iter;
		g.associated_elements(faces, vrt);

		if(faces.size() != 3)
			continue;

	//	get the connected vertices. If there are more or less than 3 connected
	//	vertices, then the connected triangles do not form a valid closed surface patch
		size_t numConVrts = 0;
		for(size_t i = 0; i < faces.size(); ++i){
			Face* f = faces[i];
			if(f->num_vertices() != 3){
				numConVrts = 0;
				break;
			}

			int baseInd = GetVertexIndex(f, vrt);

			for(size_t j = 1; j < 3; ++j){
				Vertex* cv = f->vertex((baseInd + j) % 3);
				if(cv != vrt){
					bool gotOne = false;
					for(size_t k = 0; k < numConVrts; ++k){
						if(conVrts[k] == cv){
							gotOne = true;
							break;
						}
					}
					if(!gotOne){
						if(numConVrts < 3)
							conVrts[numConVrts] = cv;
						++numConVrts;
					}
				}
			}
		}

		if(numConVrts == 3){
			if(Valence3VertexIsObsolete(aaPos[vrt], aaPos[conVrts[0]], 
										aaPos[conVrts[1]], aaPos[conVrts[2]],
										maxSquaredHeightToBaseAreaRatio))
			{
				g.create<Triangle>(TriangleDescriptor(conVrts[0], conVrts[1], conVrts[2]), faces[0]);
				g.erase(vrt);
			}
		}
	}
}


template <class TVrtIter, class TAAPos>
void ReplaceLowValenceVertices(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
							   number maxSquaredHeightToBaseAreaRatio, TAAPos aaPos)
{
	using namespace std;

	Grid::face_traits::secure_container	faces;
	Vertex* conVrts[4];

//	used for valence-4 only
	vector<Vertex*> poly;
	vector<Edge*> conEdges;

//	when we remove a vertex, the valencies of adjacent vertices change.
//	We thus operate on two candidate arrays. One which holds all vertices
//	which are current candidates, and one which records all vertices, which have
//	not been removed. After each run over alll candidates we compare both arrays.
//	If they have the same length, we're done. If not, we'll run again over all
//	new candidates.
	std::vector<Vertex*> candidates;
	std::vector<Vertex*> newCandidates;

	for(TVrtIter iter = vrtsBegin; iter != vrtsEnd; ++iter)
		newCandidates.push_back(*iter);

	while(newCandidates.size() != candidates.size())
	{
		candidates.swap(newCandidates);
		newCandidates.clear();

		for(size_t icand = 0; icand < candidates.size(); ++icand)
		{
			Vertex* vrt = candidates[icand];

		//	if vrt will be removed, we'll call 'pop_back'
			newCandidates.push_back(vrt);

			g.associated_elements(faces, vrt);

			if(faces.size() < 3 || faces.size() > 4)
				continue;

		//	get the connected vertices. If there are more or less than 3 connected
		//	vertices, then the connected triangles do not form a valid closed surface patch
			size_t numConVrts = 0;
			for(size_t i = 0; i < faces.size(); ++i){
				Face* f = faces[i];
				if(f->num_vertices() != 3){
					numConVrts = 0;
					break;
				}

				int baseInd = GetVertexIndex(f, vrt);

				for(size_t j = 1; j < 3; ++j){
					Vertex* cv = f->vertex((baseInd + j) % 3);
					if(cv != vrt){
						bool gotOne = false;
						for(size_t k = 0; k < numConVrts; ++k){
							if(conVrts[k] == cv){
								gotOne = true;
								break;
							}
						}
						if(!gotOne){
							if(numConVrts < 4)
								conVrts[numConVrts] = cv;
							++numConVrts;
						}
					}
				}
			}

			if(numConVrts != faces.size())
				continue;

			if(numConVrts == 3){
				if(Valence3VertexIsObsolete(aaPos[vrt], aaPos[conVrts[0]], 
											aaPos[conVrts[1]], aaPos[conVrts[2]],
											maxSquaredHeightToBaseAreaRatio))
				{
				//	vrt is currently the last entry in newCandidates and has to be removed
					newCandidates.pop_back();
					g.create<Triangle>(TriangleDescriptor(conVrts[0], conVrts[1], conVrts[2]), faces[0]);
					g.erase(vrt);
				}
			}
			else if(numConVrts == 4){
				conEdges.clear();
				poly.clear();
				for(size_t i = 0; i < faces.size(); ++i){
					std::pair<GridBaseObjectId, int> opObj = faces[i]->get_opposing_object(vrt);
					if(opObj.first == EDGE){
						Edge* e = g.get_edge(faces[i], opObj.second);
						if(e)
							conEdges.push_back(e);
					}
				}
				
				if(conEdges.size() != 4)
					continue;

				CreatePolyChain(poly, g, conEdges.begin(), conEdges.end());

			//	we first try the better quality-diagonal
			//	If this doesn't work, we'll try the other one.
			//	diag 0: create diagonal between vrts 0,2
			//	diag 1: between vrts 1,3
				number q[2];
				for(int diag = 0; diag < 2; ++diag){
					int i0 = diag;
					int i1 = (diag + 1);
					int i2 = (diag + 2);
					int i3 = (diag + 3) % 4;
					q[diag] = min(TriangleQuality_Area(aaPos[poly[i0]], aaPos[poly[i1]],
													   aaPos[poly[i2]]),
								  TriangleQuality_Area(aaPos[poly[i0]], aaPos[poly[i2]],
													   aaPos[poly[i3]]));
				}

				int firstDiag = 0;
				if(q[1] > q[0])
					firstDiag = 1;

				for(int idiag = 0; idiag < 2; ++idiag){
					int diag = (firstDiag + idiag) % 2;
					int i0 = diag;
					int i1 = (diag + 1);
					int i2 = (diag + 2);
					int i3 = (diag + 3) % 4;

					if(!QuadDiagonalIsValid(aaPos[poly[i0]], aaPos[poly[i1]],
											aaPos[poly[i2]], aaPos[poly[i3]]))
					{
						continue;
					}

					if(Valence4VertexIsObsolete(aaPos[vrt],
												aaPos[poly[i0]], 
												aaPos[poly[i1]],
												aaPos[poly[i2]],
												aaPos[poly[i3]],
												maxSquaredHeightToBaseAreaRatio))
					{
					//	vrt is currently the last entry in newCandidates and has to be removed
						newCandidates.pop_back();

						// FaceDescriptor fd(vrt, poly[i0], poly[i1]);
						// Face* f1 = g.get_face(fd);
						// if(!f1){
						// 	UG_LOG("num vertices in fd: " << fd.num_vertices() << endl);
						// 	UG_LOG("Couldn't find face with corners: "
						// 			<< aaPos[fd.vertex(0)]
						// 			<< ", " << aaPos[fd.vertex(1)]
						// 			<< ", " << aaPos[fd.vertex(2)] << endl);

						// 	UG_LOG("Existing faces:\n");
						// 	for(FaceIterator iter = g.begin<Face>(); iter != g.end<Face>(); ++iter){
						// 		UG_LOG(ElementDebugInfo(g, *iter) << endl);
						// 	}
						// }

						FaceDescriptor fd0(vrt, poly[i0], poly[i1]);
						FaceDescriptor fd1(vrt, poly[i2], poly[i3]);

						g.create<Triangle>(TriangleDescriptor(poly[i0],
															  poly[i1],
															  poly[i2]),
											g.get_face(fd0));
						g.create<Triangle>(TriangleDescriptor(poly[i0],
															  poly[i2],
															  poly[i3]),
											g.get_face(fd1));
						g.erase(vrt);
						break;
					}
				}
			}
		}
	}
}



// template <class TAAPosVRT, class TAANormVRT, class TAAIntVRT>
// Vertex* TryFlatRegionEdgeCollapse(Grid& grid, Edge* e,
// 							  TAAPosVRT& aaPos, TAANormVRT& aaNorm, 
// 							  TAAIntVRT& aaInt, SubsetHandler* pshMarks = NULL,
// 							  EdgeSelector* pCandidates = NULL)
// {
// 	if(pshMarks)
// 	{
// 		SubsetHandler& shMarks = *pshMarks;
// 	//	collapses are not allowed for fixed edges
// 		if(shMarks.get_subset_index(e) == REM_FIXED)
// 			return NULL;
			
// 	//	if both endpoints of are fixed vertices then
// 	//	we may not collapse
// 		int vrtSI[2];
// 		vrtSI[0] = shMarks.get_subset_index(e->vertex(0));
// 		vrtSI[1] = shMarks.get_subset_index(e->vertex(1));

// 		if((vrtSI[0] == REM_FIXED) && (vrtSI[1] == REM_FIXED))
// 			return NULL;

// 	//	if both endpoints are somehow marked, e has to be a
// 	//	crease edge
// 		if((vrtSI[0] != REM_NONE) && (vrtSI[1] != REM_NONE)
// 			&&	(shMarks.get_subset_index(e) != REM_CREASE))
// 			return NULL;
// 	}

// //	check whether the edge can be collapsed
// 	if(!EdgeCollapseIsValid(grid, e))
// 		return NULL;


// 	vector3 edgeNormal;
// 	int numNbrFaces = CalculateNormal(edgeNormal, grid, e, aaPos);
// //todo: If the edge is a crease edge, we'll have to check the crease-bending...
// 	//UG_COND_THROW(numNbrFaces > 2, "DEBUG-THROW: CREASES ARE TEMPORARILY NOT SUPPORTED IN COLLAPSE EDGE");
// 	if(numNbrFaces > 2){
// 		UG_LOG("DEBUG: IGNORING NON-MANIFOLD EDGE!\n");
// 		return NULL;
// 	}

// //	compare the normal of the edge that shall be collapsed with the normals
// //	of its corners.
// 	number cornerDotThreshold = cos(deg_to_rad<number>(0.1));

// 	number cornerDots[2];
// 	int maxCornerDotInd = 0;
// 	for(size_t i = 0; i < 2; ++i)
// 		cornerDots[i] = VecDot(edgeNormal, aaNorm[e->vertex(i)]);
// 	if(cornerDots[1] > cornerDots[0])
// 		maxCornerDotInd = 1;

// 	if(cornerDots[maxCornerDotInd] < cornerDotThreshold)
// 		return NULL;

	
// //	collapse the edge
// //todo:	When creases are involved,we have to do a better check which vertex to preserve
// 	Vertex* vrt = e->vertex(maxCornerDotInd);
// 	CollapseEdge(grid, e, vrt);

// 	return vrt;
// }


// template <class TVrtIter, class TAAPos>
// void CollapseEdgesInFlatRegions(Grid& g, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
// 							    number maxNormalDeviation, TAAPos aaPos)
// {

// }


}//	end of namespace

#endif	//__H__simplification

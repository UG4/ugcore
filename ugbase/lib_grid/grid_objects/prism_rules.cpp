/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include <cassert>
#include "prism_rules.h"
#include "rule_util.h"
#include "grid_object_ids.h"

namespace ug{
namespace prism_rules{

///	Output are the vertices of a prism rotated around its vertical axis
void RotatePrism(int vrtsOut[NUM_VERTICES], int steps)
{
	if(steps < 0)
		steps = (4 - ((-steps) % 4));

	for(int i = 0; i < 3; ++i)
		vrtsOut[(i + steps) % 3] = i;

	for(int i = 3; i < 6; ++i)
		vrtsOut[3 + (i + steps) % 3] = i;
}

//	fips the prism upside down and leaves face0 where it is
void FlipPrism(int vrtsOut[NUM_VERTICES], const int vrts[NUM_VERTICES])
{
	vrtsOut[0] = vrts[5];
	vrtsOut[1] = vrts[4];
	vrtsOut[2] = vrts[3];
	vrtsOut[3] = vrts[2];
	vrtsOut[4] = vrts[1];
	vrtsOut[5] = vrts[0];
}


int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut, vector3*, bool* isSnapPoint)
{
	newCenterOut = false;
//	If a refinement rule is not implemented, fillCount will stay at 0.
//	Before returning, we'll check, whether fillCount is at 0 and will
//	perform refinement with an additional vertex in this case.

//	corner status is used to mark vertices, which are connected to refined edges
//	the value tells, how many associated edges are refined.
	int cornerStatus[6] = {0, 0, 0, 0, 0, 0};

//	here we'll store the index of each edge, which will be refined.
//	Size of the valid sub-array is numNewVrts, which is calculated below.
	int refEdgeInds[NUM_EDGES];

//	count the number of new vertices and fill newVrtEdgeInds
	int numNewVrts = 0;
	for(int i = 0; i < NUM_EDGES; ++i){
		if(newEdgeVrts[i]){
			refEdgeInds[numNewVrts] = i;
			++numNewVrts;

		//	adjust corner status of associated vertices
			const int* evi = EDGE_VRT_INDS[i];
			++cornerStatus[evi[0]];
			++cornerStatus[evi[1]];
		}
	}

//	snap-point handling
	int numSnapPoints = 0;
	int snapPointIndex[NUM_VERTICES];
	if(isSnapPoint){
		for(int i = 0; i < NUM_VERTICES; ++i){
			if(isSnapPoint[i]){
				snapPointIndex[numSnapPoints] = i;
				++numSnapPoints;
			}
		}
	}

//	set to true if the algorithm already tried to find a special rule for the
//	specified snap-points. Always true if no snap points were provided.
	bool snapPointsProcessed = (numSnapPoints == 0);

//	the fillCount tells how much data has already been written to newIndsOut.
	int fillCount = 0;

//	convenience - indices where new edge-vrts, new face-vrts and new vol-vrts begin.
	constexpr int E = NUM_VERTICES;
	constexpr int F = NUM_VERTICES + NUM_EDGES;
//	const int V = NUM_VERTICES + NUM_EDGES + NUM_FACES;

//	depending on the number of new vertices, we will now apply different
//	refinement rules. Further distinction may have to be done.
	switch(numNewVrts){
		case 0:
		{
		//	simply put the default prism back to newIndsOut
			newIndsOut[fillCount++] = GridObjectID::GOID_PRISM;
			newIndsOut[fillCount++] = 0;
			newIndsOut[fillCount++] = 1;
			newIndsOut[fillCount++] = 2;
			newIndsOut[fillCount++] = 3;
			newIndsOut[fillCount++] = 4;
			newIndsOut[fillCount++] = 5;
		}break;

		case 1:
		{
		//	we'll iterate over all faces. If a face does not contain the
		//	refined edge, then we'll create a tetrahedron or pyramid from it.
			const int refEdge = refEdgeInds[0];
			const int nVrt = refEdge + E;

			int& fi = fillCount;
			int* inds = newIndsOut;

			if(numSnapPoints == 2){
				const int snapEdge = EDGE_FROM_VRTS[snapPointIndex[0]][snapPointIndex[1]];
				if(snapEdge != -1
					&& !IS_SIDE_EDGE[snapEdge]
					&& IS_SIDE_EDGE[refEdge]
					&& !isSnapPoint[EDGE_VRT_INDS[refEdge][0]]
					&& !isSnapPoint[EDGE_VRT_INDS[refEdge][1]])
				{
				//	find the quad which contains the snapEdge
					int snapQuad = -1;
					for(int iquad = 0; iquad < 3; ++iquad){
						if(FACE_CONTAINS_EDGE[QUADS[iquad]][snapEdge]){
							snapQuad = QUADS[iquad];
							break;
						}
					}

					UG_ASSERT(snapQuad != -1, "The corresponding quad has to exist!");

					int p[NUM_VERTICES];
					if(IS_BOTTOM_EDGE[snapEdge]){
						int tp[NUM_VERTICES];
						RotatePrism(tp, 3 - snapQuad);
						FlipPrism(p, tp);
					}
					else
						RotatePrism(p, 3 - snapQuad);


					// UG_LOG("snapEdge: " << snapEdge << "\n");
					// UG_LOG("refEdge: " << refEdge << "\n");
					// UG_LOG("snapPointIndex[0]: " << snapPointIndex[0] << "\n");
					// UG_LOG("snapPointIndex[1]: " << snapPointIndex[1] << "\n");
					// UG_LOG("snapQuad: " << snapQuad << "\n");
					
					// UG_LOG("p: ");
					// for(int ivrt = 0; ivrt < NUM_VERTICES; ++ivrt){
					// 	UG_LOG(p[ivrt] << ", ");
					// }
					// UG_LOG(std::endl);



					const int e14 = E + EDGE_FROM_VRTS[p[1]][p[4]];

					//	we have to create two new prisms
					int& fi = fillCount;
					int* inds = newIndsOut;

					inds[fi++] = GridObjectID::GOID_PRISM;
					inds[fi++] = p[0];	inds[fi++] = p[1];	inds[fi++] = p[2];
					inds[fi++] = p[3];	inds[fi++] = e14;	inds[fi++] = p[5];

					inds[fi++] = GridObjectID::GOID_TETRAHEDRON;
					inds[fi++] = p[3];	inds[fi++] = e14;	inds[fi++] = p[5];
					inds[fi++] = p[4];

					snapPointsProcessed = true;
				}
			}

			if(numSnapPoints == 0 || !snapPointsProcessed){
				for(int i = 0; i < NUM_FACES; ++i){
					if(!FACE_CONTAINS_EDGE[i][refEdge]){
						const int* f = FACE_VRT_INDS[i];

					//	if f is a tri, then we'll create a tetrahedron, if it is a
					//	quad, then we'll create a pyramid.
						if(f[3] == -1){
						//	create a tetrahedron
							inds[fi++] = GridObjectID::GOID_TETRAHEDRON;
							inds[fi++] = f[0];	inds[fi++] = f[1];
							inds[fi++] = f[2];	inds[fi++] = nVrt;
						}
						else{
						//	create a prism
							inds[fi++] = GridObjectID::GOID_PYRAMID;
							inds[fi++] = f[0];	inds[fi++] = f[1];
							inds[fi++] = f[2];	inds[fi++] = f[3];
							inds[fi++] = nVrt;
						}
					}
				}
			}
		}break;

		case 2:
		{
		//	three important cases exist
			const int re0 = refEdgeInds[0];
			const int re1 = refEdgeInds[1];
		//	index of the face that contains both edges. This may be -1.
			const int refFaceInd = FACE_FROM_EDGES[re0][re1];

		//	we currently don't handle cases where refFaceInd == -1
			if(refFaceInd == -1)
				break;

			if(FACE_VRT_INDS[refFaceInd][3] != -1){
			//	a side-face (quad) is refined
			//	rotate the prism, so that the refined face lies on face 3.
				int p[NUM_VERTICES];
				RotatePrism(p, 3 - refFaceInd);

				int& fi = fillCount;
				int* inds = newIndsOut;

				const int e3 = EDGE_FROM_VRTS[p[0]][p[3]] + E;
				const int e5 = EDGE_FROM_VRTS[p[2]][p[5]] + E;

			//	check whether the cut is horizontal or vertical (ignore others)
				if(IS_SIDE_EDGE[re0] && IS_SIDE_EDGE[re1]){
				//	horizontal cut
					if(numSnapPoints == 1){
					//	build a pyramid and a prism
						if(isSnapPoint[p[1]]){
							inds[fi++] = GridObjectID::GOID_PYRAMID;
							inds[fi++] = p[0];	inds[fi++] = p[2];	inds[fi++] = e5;
							inds[fi++] = e3;	inds[fi++] = p[1];

							inds[fi++] = GridObjectID::GOID_PRISM;
							inds[fi++] = e3;	inds[fi++] = p[1];	inds[fi++] = e5;
							inds[fi++] = p[3];	inds[fi++] = p[4];	inds[fi++] = p[5];
							snapPointsProcessed = true;
						}
						else if(isSnapPoint[p[4]]){
							inds[fi++] = GridObjectID::GOID_PYRAMID;
							inds[fi++] = e3;	inds[fi++] = e5;	inds[fi++] = p[5];
							inds[fi++] = p[3];	inds[fi++] = p[4];

							inds[fi++] = GridObjectID::GOID_PRISM;
							inds[fi++] = p[0];	inds[fi++] = p[1];	inds[fi++] = p[2];
							inds[fi++] = e3;	inds[fi++] = p[4];	inds[fi++] = e5;
							snapPointsProcessed = true;
						}
					}

					if(numSnapPoints == 0 || !snapPointsProcessed){
					//	build two pyramids and a tetrahedron.
						inds[fi++] = GridObjectID::GOID_PYRAMID;
						inds[fi++] = p[0];	inds[fi++] = p[2];	inds[fi++] = e5;
						inds[fi++] = e3;	inds[fi++] = p[1];

						inds[fi++] = GridObjectID::GOID_PYRAMID;
						inds[fi++] = e3;	inds[fi++] = e5;	inds[fi++] = p[5];
						inds[fi++] = p[3];	inds[fi++] = p[4];

						inds[fi++] = GridObjectID::GOID_TETRAHEDRON;
						inds[fi++] = p[1];	inds[fi++] = e3;
						inds[fi++] = p[4];	inds[fi++] = e5;
					}
				}
				else if(!IS_SIDE_EDGE[re0] && !IS_SIDE_EDGE[re1]){
				//	vertical cut
				//	create two prisms
					const int e2 = EDGE_FROM_VRTS[p[0]][p[2]] + E;
					const int e8 = EDGE_FROM_VRTS[p[3]][p[5]] + E;

					inds[fi++] = GridObjectID::GOID_PRISM;
					inds[fi++] = p[0];	inds[fi++] = p[1];	inds[fi++] = e2;
					inds[fi++] = p[3];	inds[fi++] = p[4];	inds[fi++] = e8;
					inds[fi++] = GridObjectID::GOID_PRISM;
					inds[fi++] = p[1];	inds[fi++] = p[2];	inds[fi++] = e2;
					inds[fi++] = p[4];	inds[fi++] = p[5];	inds[fi++] = e8;
				}
			}
			else{
			//	refined face is a triangle - either top or bottom.
			//	some cases still missing...
				if(numSnapPoints == 2){
				//	get the non-marked edge of refFace
					int nonRefEdge = -1;
					for(int iedge = 0; iedge < 3; ++iedge){
						if(!newEdgeVrts[FACE_EDGE_INDS[refFaceInd][iedge]]){
							nonRefEdge = FACE_EDGE_INDS[refFaceInd][iedge];
							break;
						}
					}
					UG_ASSERT(nonRefEdge == -1, "A non-ref edge has to exist!");

				//	get the edge which connects the two snap points
					int snapEdge = EDGE_FROM_VRTS[snapPointIndex[0]][snapPointIndex[1]];
					if(snapEdge == -1 || IS_SIDE_EDGE[snapEdge])
						break;

					int snapQuad = FACE_FROM_EDGES[nonRefEdge][snapEdge];
					if(snapQuad == -1 || snapQuad == refFaceInd)
						break;


					int p[NUM_VERTICES];
					if(refFaceInd == BOTTOM_FACE){
						int tp[NUM_VERTICES];
						RotatePrism(tp, 3 - snapQuad);
						FlipPrism(p, tp);
					}
					else
						RotatePrism(p, 3 - snapQuad);

					const int e34 = E + EDGE_FROM_VRTS[p[3]][p[4]];
					const int e45 = E + EDGE_FROM_VRTS[p[4]][p[5]];
					//	we have to create two new prisms
					int& fi = fillCount;
					int* inds = newIndsOut;

					inds[fi++] = GridObjectID::GOID_PRISM;
					inds[fi++] = p[0];	inds[fi++] = p[3];	inds[fi++] = e34;
					inds[fi++] = p[2];	inds[fi++] = p[5];	inds[fi++] = e45;

					inds[fi++] = GridObjectID::GOID_PRISM;
					inds[fi++] = p[0];	inds[fi++] = p[1];	inds[fi++] = p[2];
					inds[fi++] = e34;	inds[fi++] = p[4];	inds[fi++] = e45;

					snapPointsProcessed = true;
				}
			}
		}break;

		case 3:
		{
		//	here we currently only handle a horizontal slice
			if(IS_SIDE_EDGE[refEdgeInds[0]] && IS_SIDE_EDGE[refEdgeInds[1]]
			   && IS_SIDE_EDGE[refEdgeInds[2]])
			{
				constexpr int e3 = E + 3;
				constexpr int e4 = E + 4;
				constexpr int e5 = E + 5;

			//	we have to create two new prisms
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = 0;		inds[fi++] = 1;		inds[fi++] = 2;
				inds[fi++] = e3;	inds[fi++] = e4;	inds[fi++] = e5;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = e3;	inds[fi++] = e4;	inds[fi++] = e5;
				inds[fi++] = 3;		inds[fi++] = 4;		inds[fi++] = 5;
			}
		}break;

		case 4:
		{
		//	only a nice, vertical cut is directly supported in the moment
		//	this means that 2 marked edges are top-edges and two are bottom edges
			int numTop = 0;
			int numBottom = 0;
			for(int i = 0; i < 4; ++i){
				if(IS_TOP_EDGE[refEdgeInds[i]])
					++numTop;
				else if(IS_BOTTOM_EDGE[refEdgeInds[i]])
					++numBottom;
			}

			if((numTop == 2) && (numBottom == 2)){
			//	find a quadrilateral, whose corner states are all 1.
				int freeFaceInd = -1;
				for(int i = 0; i < NUM_QUADS; ++i){
					const int* f = FACE_VRT_INDS[QUADS[i]];
					bool allOne = true;
					for(int j = 0; j < 4; ++j){
						if(cornerStatus[f[j]] != 1){
							allOne = false;
							break;
						}
					}

					if(allOne){
						freeFaceInd = QUADS[i];
						break;
					}
				}

				if(freeFaceInd != -1){
				//	we got the free quad. It is thus clear, that all associated
				//	edges, which do not lie in the quad itself, will be refined.
				//	create a hexahedron and a prism.
				//	First however rotate the prism so that the free quad is at
				//	index 3.
					int p[NUM_VERTICES];
					RotatePrism(p, 3 - freeFaceInd);

				//	some important indices.
					const int e0 = EDGE_FROM_VRTS[p[0]][p[1]] + E;
					const int e1 = EDGE_FROM_VRTS[p[1]][p[2]] + E;
					const int e6 = EDGE_FROM_VRTS[p[3]][p[4]] + E;
					const int e7 = EDGE_FROM_VRTS[p[4]][p[5]] + E;

				//	now create the elements
					int& fi = fillCount;
					int* inds = newIndsOut;

					inds[fi++] = GridObjectID::GOID_HEXAHEDRON;
					inds[fi++] = p[0];	inds[fi++] = e0;
					inds[fi++] = e1;	inds[fi++] = p[2];
					inds[fi++] = p[3];	inds[fi++] = e6;
					inds[fi++] = e7;	inds[fi++] = p[5];

					inds[fi++] = GridObjectID::GOID_PRISM;
					inds[fi++] = e0;	inds[fi++] = p[1];	inds[fi++] = e1;
					inds[fi++] = e6;	inds[fi++] = p[4];	inds[fi++] = e7;
				}
			}
		}break;

		case 6:
		{
		//	we currently only support the case, where only top and bottom edges
		//	are refined, resulting in 4 new prisms.
			bool sideRefined = false;
			for(int i = 0; i < 6; ++i){
				if(IS_SIDE_EDGE[refEdgeInds[i]]){
					sideRefined = true;
					break;
				}
			}

			if(!sideRefined){
			//	create 4 new prisms. First however store the new edge vertices.
			constexpr int e0 = E;
			constexpr int e1 = E + 1;
			constexpr int e2 = E + 2;
			constexpr int e6 = E + 6;
			constexpr int e7 = E + 7;
			constexpr int e8 = E + 8;

				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = 0;		inds[fi++] = e0;	inds[fi++] = e2;
				inds[fi++] = 3;		inds[fi++] = e6;	inds[fi++] = e8;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = 1;		inds[fi++] = e1;	inds[fi++] = e0;
				inds[fi++] = 4;		inds[fi++] = e7;	inds[fi++] = e6;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = 2;		inds[fi++] = e2;	inds[fi++] = e1;
				inds[fi++] = 5;		inds[fi++] = e8;	inds[fi++] = e7;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = e0;	inds[fi++] = e1;	inds[fi++] = e2;
				inds[fi++] = e6;	inds[fi++] = e7;	inds[fi++] = e8;
			}
		}break;

		case 7:
		{
		//	only the case in which two quads are completly refined is regarded.
		//	In this case corresponding edges of the bottom and top triangle are
		//	marked/unmarked.
		//	order bv so, that the one connected to two marked base-edges
		//	is bv[0].
			int bv[3], tv[3];
			bool twoQuadsMarked = false;
			for(int iedge = 0; iedge < 3; ++iedge){
				if(!newEdgeVrts[iedge]){
					if(!newEdgeVrts[iedge + 6]){
						twoQuadsMarked = true;

						bv[0] = (iedge + 2) % 3;
						bv[1] = (iedge + 3) % 3;
						bv[2] = (iedge + 4) % 3;

						tv[0] = bv[0] + 3;
						tv[1] = bv[1] + 3;
						tv[2] = bv[2] + 3;
					}
					break;
				}
			}


			if(twoQuadsMarked){
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = bv[0];		inds[fi++] = E + bv[0];
				inds[fi++] = E + bv[2];
				
				inds[fi++] = E + bv[0] + 3;		inds[fi++] = F + QUADS[bv[0]];
				inds[fi++] = F + QUADS[bv[2]];


				inds[fi++] = GridObjectID::GOID_HEXAHEDRON;
				inds[fi++] = E + bv[0]; inds[fi++] = bv[1];
				inds[fi++] = bv[2];		inds[fi++] = E + bv[2];

				inds[fi++] = F + QUADS[bv[0]];	inds[fi++] = E + bv[1] + 3;
				inds[fi++] = E + bv[2] + 3;		inds[fi++] = F + QUADS[bv[2]];

				inds[fi++] = GridObjectID::GOID_PRISM;
				inds[fi++] = E + bv[0] + 3;	inds[fi++] = F + QUADS[bv[0]];
				inds[fi++] = F + QUADS[bv[2]];

				inds[fi++] = tv[0];		inds[fi++] = E + 3 + tv[0];
				inds[fi++] = E + 3 + tv[2];


				inds[fi++] = GridObjectID::GOID_HEXAHEDRON;
				inds[fi++] = F + QUADS[bv[0]];	inds[fi++] = E + bv[1] + 3;
				inds[fi++] = E + bv[2] + 3;	inds[fi++] = F + QUADS[bv[2]];

				inds[fi++] = E + 3 + tv[0]; inds[fi++] = tv[1];
				inds[fi++] = tv[2];			inds[fi++] = E + 3 + tv[2];
			}
		}break;

		case 9:	// regular refine
		{
		//	we have to create 8 new prisms
			int& fi = fillCount;
			int* inds = newIndsOut;
			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = 0;			inds[fi++] = E;			inds[fi++] = E + 2;
			inds[fi++] = E + 3;		inds[fi++] = F + 1;		inds[fi++] = F + 3;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = E;			inds[fi++] = 1;			inds[fi++] = E + 1;
			inds[fi++] = F + 1;		inds[fi++] = E + 4;		inds[fi++] = F + 2;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = E + 2;		inds[fi++] = E + 1;		inds[fi++] = 2;
			inds[fi++] = F + 3;		inds[fi++] = F + 2;		inds[fi++] = E + 5;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = E + 1;		inds[fi++] = E + 2;		inds[fi++] = E;
			inds[fi++] = F + 2;		inds[fi++] = F + 3;		inds[fi++] = F + 1;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = E + 3;		inds[fi++] = F + 1;		inds[fi++] = F + 3;
			inds[fi++] = 3;			inds[fi++] = E + 6;		inds[fi++] = E + 8;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = F + 1;		inds[fi++] = E + 4;		inds[fi++] = F + 2;
			inds[fi++] = E + 6;		inds[fi++] = 4;			inds[fi++] = E + 7;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = F + 3;		inds[fi++] = F + 2;		inds[fi++] = E + 5;
			inds[fi++] = E + 8;		inds[fi++] = E + 7;		inds[fi++] = 5;

			inds[fi++] = GridObjectID::GOID_PRISM;
			inds[fi++] = F + 2;		inds[fi++] = F + 3;		inds[fi++] = F + 1;
			inds[fi++] = E + 7;		inds[fi++] = E + 8;		inds[fi++] = E + 6;
		}break;
	}

	if(!snapPointsProcessed){
		UG_LOG("WARNING: Invalid or unsupported snap-point distribution detected. Ignoring snap-points for this element.\n");
	}

	if(fillCount == 0){
	//	call recursive refine
		fillCount = shared_rules::RecursiveRefine(newIndsOut, newEdgeVrts,
									FACE_VRT_INDS, FACE_EDGE_INDS, NUM_VERTICES,
									NUM_EDGES, NUM_FACES);

	//	the rule requires a new center vertex
		newCenterOut = true;
	}

	return fillCount;
}


int CollapseEdge (int* newIndsOut, int v0, int v1)
{
//	get the edge which is connected by v0 and v1
	const int edgeInd = EDGE_FROM_VRTS[v0][v1];

//	if the edge doesn't exist or if one of the triangle-edges is collapsed,
//	no volume can be created
	if((edgeInd < 3) || (edgeInd > 5))
		return 0;

//	the index of the opposing edge of the base quadrilateral
//	the two non-collapsed vertical edges
	const int opEdgeInd0 = 3 + (edgeInd + 1) % 3;
	const int opEdgeInd1 = 3 + (edgeInd + 2) % 3;

//	the base of the pyramid
	const int baseFaceInd = FACE_FROM_EDGES[opEdgeInd0][opEdgeInd1];

//	the resulting volume is a pyramid
	int fi = 0;
	newIndsOut[fi++] = GridObjectID::GOID_PYRAMID;
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][0];
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][1];
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][2];
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][3];
	newIndsOut[fi++] = v0;

	return fi;
}

bool IsRegularRefRule(const int edgeMarks)
{
	static constexpr int allEdges = 511;	// 0b111111111
	static constexpr int allHEdges = 455;	// 0b111000111
	static constexpr int allVEdges = 56;	// 0b000111000

	return		(edgeMarks == allEdges)
			||	(edgeMarks == allHEdges)
			||	(edgeMarks == allVEdges);
}

}// end of namespace
}// end of namespace

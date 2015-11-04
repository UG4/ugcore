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


int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut, vector3*)
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

//	the fillCount tells how much data has already been written to newIndsOut.
	int fillCount = 0;

//	convenience - indices where new edge-vrts, new face-vrts and new vol-vrts begin.
	const int E = NUM_VERTICES;
	const int F = NUM_VERTICES + NUM_EDGES;
//	const int V = NUM_VERTICES + NUM_EDGES + NUM_FACES;

//	depending on the number of new vertices, we will now apply different
//	refinement rules. Further distinction may have to be done.
	switch(numNewVrts){
		case 0:
		{
		//	simply put the default prism back to newIndsOut
			newIndsOut[fillCount++] = GOID_PRISM;
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

			for(int i = 0; i < NUM_FACES; ++i){
				if(!FACE_CONTAINS_EDGE[i][refEdge]){
					const int* f = FACE_VRT_INDS[i];

				//	if f is a tri, then we'll create a tetrahedron, if it is a
				//	quad, then we'll create a pyramid.
					if(f[3] == -1){
					//	create a tetrahedron
						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = f[0];	inds[fi++] = f[1];
						inds[fi++] = f[2];	inds[fi++] = nVrt;
					}
					else{
					//	create a prism
						inds[fi++] = GOID_PYRAMID;
						inds[fi++] = f[0];	inds[fi++] = f[1];
						inds[fi++] = f[2];	inds[fi++] = f[3];
						inds[fi++] = nVrt;
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

			//	check whether the cut is horizontal or vertical (ignore others)
				if(IS_SIDE_EDGE[re0] && IS_SIDE_EDGE[re1]){
				//	horizontal cut
				//	build two pyramids and a tetrahedron.
					const int e3 = EDGE_FROM_VRTS[p[0]][p[3]] + E;
					const int e5 = EDGE_FROM_VRTS[p[2]][p[5]] + E;

					inds[fi++] = GOID_PYRAMID;
					inds[fi++] = p[0];	inds[fi++] = p[2];	inds[fi++] = e5;
					inds[fi++] = e3;	inds[fi++] = p[1];

					inds[fi++] = GOID_PYRAMID;
					inds[fi++] = e3;	inds[fi++] = e5;	inds[fi++] = p[5];
					inds[fi++] = p[3];	inds[fi++] = p[4];

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = p[1];	inds[fi++] = e3;
					inds[fi++] = p[4];	inds[fi++] = e5;
				}
				else if(!IS_SIDE_EDGE[re0] && !IS_SIDE_EDGE[re1]){
				//	vertical cut
				//	create two prisms
					const int e2 = EDGE_FROM_VRTS[p[0]][p[2]] + E;
					const int e8 = EDGE_FROM_VRTS[p[3]][p[5]] + E;

					inds[fi++] = GOID_PRISM;
					inds[fi++] = p[0];	inds[fi++] = p[1];	inds[fi++] = e2;
					inds[fi++] = p[3];	inds[fi++] = p[4];	inds[fi++] = e8;
					inds[fi++] = GOID_PRISM;
					inds[fi++] = p[1];	inds[fi++] = p[2];	inds[fi++] = e2;
					inds[fi++] = p[4];	inds[fi++] = p[5];	inds[fi++] = e8;
				}
			}
			else{
			//	refined face is a triangle - either top or bottom.
			//	this should be implemented somewhen...
			}
		}break;

		case 3:
		{
		//	here we currently only handle a horizontal slice
			if(IS_SIDE_EDGE[refEdgeInds[0]] && IS_SIDE_EDGE[refEdgeInds[1]]
			   && IS_SIDE_EDGE[refEdgeInds[2]])
			{
				const int e3 = E + 3;
				const int e4 = E + 4;
				const int e5 = E + 5;

			//	we have to create two new prisms
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = 0;		inds[fi++] = 1;		inds[fi++] = 2;
				inds[fi++] = e3;	inds[fi++] = e4;	inds[fi++] = e5;

				inds[fi++] = GOID_PRISM;
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

					inds[fi++] = GOID_HEXAHEDRON;
					inds[fi++] = p[0];	inds[fi++] = e0;
					inds[fi++] = e1;	inds[fi++] = p[2];
					inds[fi++] = p[3];	inds[fi++] = e6;
					inds[fi++] = e7;	inds[fi++] = p[5];

					inds[fi++] = GOID_PRISM;
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
				const int e0 = E;
				const int e1 = E + 1;
				const int e2 = E + 2;
				const int e6 = E + 6;
				const int e7 = E + 7;
				const int e8 = E + 8;

				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = 0;		inds[fi++] = e0;	inds[fi++] = e2;
				inds[fi++] = 3;		inds[fi++] = e6;	inds[fi++] = e8;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = 1;		inds[fi++] = e1;	inds[fi++] = e0;
				inds[fi++] = 4;		inds[fi++] = e7;	inds[fi++] = e6;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = 2;		inds[fi++] = e2;	inds[fi++] = e1;
				inds[fi++] = 5;		inds[fi++] = e8;	inds[fi++] = e7;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = e0;	inds[fi++] = e1;	inds[fi++] = e2;
				inds[fi++] = e6;	inds[fi++] = e7;	inds[fi++] = e8;
			}
		}break;

		case 7:
		{
		//	only the case in which two quads are completly refined is regarded.

		//	todo ...
		}break;

		case 9:	// regular refine
		{
		//	we have to create 8 new prisms
			int& fi = fillCount;
			int* inds = newIndsOut;
			inds[fi++] = GOID_PRISM;
			inds[fi++] = 0;			inds[fi++] = E;			inds[fi++] = E + 2;
			inds[fi++] = E + 3;		inds[fi++] = F + 1;		inds[fi++] = F + 3;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = E;			inds[fi++] = 1;			inds[fi++] = E + 1;
			inds[fi++] = F + 1;		inds[fi++] = E + 4;		inds[fi++] = F + 2;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = E + 2;		inds[fi++] = E + 1;		inds[fi++] = 2;
			inds[fi++] = F + 3;		inds[fi++] = F + 2;		inds[fi++] = E + 5;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = E + 1;		inds[fi++] = E + 2;		inds[fi++] = E;
			inds[fi++] = F + 2;		inds[fi++] = F + 3;		inds[fi++] = F + 1;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = E + 3;		inds[fi++] = F + 1;		inds[fi++] = F + 3;
			inds[fi++] = 3;			inds[fi++] = E + 6;		inds[fi++] = E + 8;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = F + 1;		inds[fi++] = E + 4;		inds[fi++] = F + 2;
			inds[fi++] = E + 6;		inds[fi++] = 4;			inds[fi++] = E + 7;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = F + 3;		inds[fi++] = F + 2;		inds[fi++] = E + 5;
			inds[fi++] = E + 8;		inds[fi++] = E + 7;		inds[fi++] = 5;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = F + 2;		inds[fi++] = F + 3;		inds[fi++] = F + 1;
			inds[fi++] = E + 7;		inds[fi++] = E + 8;		inds[fi++] = E + 6;
		}break;
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
	newIndsOut[fi++] = GOID_PYRAMID;
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][0];
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][1];
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][2];
	newIndsOut[fi++] = FACE_VRT_INDS[baseFaceInd][3];
	newIndsOut[fi++] = v0;

	return fi;
}

}// end of namespace
}// end of namespace

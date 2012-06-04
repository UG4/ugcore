// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.05.2011 (m,d,y)

#include <cassert>
#include "pyramid_rules.h"
#include "rule_util.h"

namespace ug{
namespace pyra_rules{

///	Output are the vertices of a pyramid rotated around its vertical axis
void RotatePyramid(int vrtsOut[NUM_VERTICES], int steps)
{
	if(steps < 0)
		steps = (4 - ((-steps) % 4));

	for(int i = 0; i < 4; ++i)
		vrtsOut[(i + steps) % 4] = i;

	vrtsOut[4] = 4;
}




int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut, vector3*)
{
	newCenterOut = false;
//	If a refinement rule is not implemented, fillCount will stay at 0.
//	Before returning, we'll check, whether fillCount is at 0 and will
//	perform refinement with an additional vertex in this case.

//	corner status is used to mark vertices, which are connected to refined edges
//	the value tells, how many associated edges are refined.
	int cornerStatus[5] = {0, 0, 0, 0, 0};

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
		//	simply put the default pyramid back to newIndsOut
			newIndsOut[fillCount++] = 5;
			newIndsOut[fillCount++] = 0;
			newIndsOut[fillCount++] = 1;
			newIndsOut[fillCount++] = 2;
			newIndsOut[fillCount++] = 3;
			newIndsOut[fillCount++] = 4;
		}break;

		case 1: // only one edge has to be refined.
		{
			const int refEdge = refEdgeInds[0];
			const int nVrt = refEdge + E;

			int& fi = fillCount;
			int* inds = newIndsOut;

		//	iterate through the faces. If a face does not contain refEdge,
		//	then we'll create a tetrahedron or a pyramid from the face to the
		//	new vertex
			for(int i = 0; i < NUM_FACES; ++i){
				if(!FACE_CONTAINS_EDGE[i][refEdge]){
				//	if the face is a triangle, then create a tetrahedron.
				//	if it is a quadrilateral, then create a pyramid.
					const int* f = FACE_VRT_INDS[i];

					if(f[3] == -1){
					//	create a tet
						inds[fi++] = 4;
						inds[fi++] = f[0];	inds[fi++] = f[1];
						inds[fi++] = f[2];	inds[fi++] = nVrt;
					}
					else{
					//	create a pyramid
						inds[fi++] = 5;
						inds[fi++] = f[0];	inds[fi++] = f[1];
						inds[fi++] = f[2];	inds[fi++] = f[3];
						inds[fi++] = nVrt;
					}
				}
			}
		}break;

		case 2: // multiple cases exist. Only some are directly supported.
		{
		//	get the face in which the refEdges lie (if there is one at all)
			int refFaceInd = FACE_FROM_EDGES[refEdgeInds[0]][refEdgeInds[1]];

			if(refFaceInd != -1){
			//	both refined edges lie in one face.
			//	they either lie in the quad or in one of the triangles
				const int* f = FACE_VRT_INDS[refFaceInd];
				if(f[3] != -1){
				//	a quad
				//	we have to check whether the vertices lie at opposing
				//	edges or at adjacent edges. We'll check corner-states for that.
				//	if there is a corner with status 2, they lie on adjacent ones.
					int corner2 = -1;
					int corner2LocInd = -1;
					for(int i = 0; i < 4; ++i){
						if(cornerStatus[f[i]] == 2){
							corner2 = f[i];
							corner2LocInd = i;
							break;
						}
					}

					if(corner2 == -1){
					//	opposing edges. Rotate the pyramid so that the edge
					//	between v0 and v1 will not be splitted.
						int steps = 0;
						int firstEdgeInd = EDGE_FROM_VRTS[0][1];
						if(newEdgeVrts[firstEdgeInd] != 0)
							steps = 1;

					//	the rotated pyramid
						int p[5];
						RotatePyramid(p, steps);

					//	important vertex indices
						int v1v2 = EDGE_FROM_VRTS[p[1]][p[2]] + E;
						int v3v0 = EDGE_FROM_VRTS[p[3]][p[0]] + E;

					//	now build two pyramids
						int& fi = fillCount;
						int* inds = newIndsOut;

						inds[fi++] = 5;
						inds[fi++] = p[0];		inds[fi++] = p[1];
						inds[fi++] = v1v2;		inds[fi++] = v3v0;
						inds[fi++] = p[4];

						inds[fi++] = 5;
						inds[fi++] = v3v0;		inds[fi++] = v1v2;
						inds[fi++] = p[2];		inds[fi++] = p[3];
						inds[fi++] = p[4];
					}
					else{
					//	adjacent edges. Rotate the pyramid so that the corner is
					//	at vertex 2
						int p[5];
						int steps = 2 - corner2LocInd;
						RotatePyramid(p, steps);

					//	important vertex indices
						int v1v2 = EDGE_FROM_VRTS[p[1]][p[2]] + E;
						int v2v3 = EDGE_FROM_VRTS[p[2]][p[3]] + E;

					//	now we have to create 4 tetrahedrons
						int& fi = fillCount;
						int* inds = newIndsOut;

						inds[fi++] = 4;
						inds[fi++] = p[0];		inds[fi++] = p[1];
						inds[fi++] = v1v2;		inds[fi++] = p[4];

						inds[fi++] = 4;
						inds[fi++] = v1v2;		inds[fi++] = p[2];
						inds[fi++] = v2v3;		inds[fi++] = p[4];

						inds[fi++] = 4;
						inds[fi++] = v2v3;		inds[fi++] = p[3];
						inds[fi++] = p[0];		inds[fi++] = p[4];

						inds[fi++] = 4;
						inds[fi++] = v1v2;		inds[fi++] = v2v3;
						inds[fi++] = p[0];		inds[fi++] = p[4];
					}
				}
				else{
				//	a tri.
				//	only nice cuts are currently directly supported
					if(cornerStatus[4] == 2){
					//	we'll create a prism and a pyramid. However, first
					//	we'll rotate the pyramid, so that the cut lies in face 2
						int p[5];
						RotatePyramid(p, (2 - refFaceInd));

						int e1 = EDGE_FROM_VRTS[p[1]][p[4]] + E;
						int e2 = EDGE_FROM_VRTS[p[2]][p[4]] + E;

						int& fi = fillCount;
						int* inds = newIndsOut;

						inds[fi++] = 6;
						inds[fi++] = p[1];	inds[fi++] = p[0];	inds[fi++] = e1;
						inds[fi++] = p[2];	inds[fi++] = p[3];	inds[fi++] = e2;

						inds[fi++] = 5;
						inds[fi++] = p[0];	inds[fi++] = e1;	inds[fi++] = e2;
						inds[fi++] = p[3];	inds[fi++] = p[4];
					}
				}
			}
		}break;

		case 3:
		{
		//	currently we only directly support refinement of the base
		//	get the face in which the refEdges lie (if there is one at all)
			int refFaceInd = FACE_FROM_EDGES[refEdgeInds[0]][refEdgeInds[1]];
			int tFaceInd = FACE_FROM_EDGES[refEdgeInds[1]][refEdgeInds[2]];

			if(refFaceInd == tFaceInd && refFaceInd != -1){
				const int* f = FACE_VRT_INDS[refFaceInd];
			//	make sure that the face is a quadrilateral
				if(f[3] != -1){
				//	rotate the pyramid, so that the edge v0v1 is not refined
					int freeEdge = -1;
					for(int i = 0; i < 4; ++i){
						int e = FACE_EDGE_INDS[refFaceInd][i];
						assert(e != -1);
						if(!newEdgeVrts[e]){
							freeEdge = e;
							break;
						}
					}

					assert(freeEdge != -1);

				//	the rotated pyramid
					int p[5];
					RotatePyramid(p, -freeEdge);

				//	important vertex indices
					int v1v2 = EDGE_FROM_VRTS[p[1]][p[2]] + E;
					int v2v3 = EDGE_FROM_VRTS[p[2]][p[3]] + E;
					int v3v0 = EDGE_FROM_VRTS[p[3]][p[0]] + E;

				//	now we have to create 1 pyramid and 3 tetrahedrons
					int& fi = fillCount;
					int* inds = newIndsOut;

					inds[fi++] = 5;
					inds[fi++] = p[0];	inds[fi++] = p[1];	inds[fi++] = v1v2;
					inds[fi++] = v3v0;	inds[fi++] = p[4];

					inds[fi++] = 4;
					inds[fi++] = v1v2;	inds[fi++] = p[2];
					inds[fi++] = v2v3;	inds[fi++] = p[4];

					inds[fi++] = 4;
					inds[fi++] = v2v3;	inds[fi++] = p[3];
					inds[fi++] = v3v0;	inds[fi++] = p[4];

					inds[fi++] = 4;
					inds[fi++] = v1v2;	inds[fi++] = v2v3;
					inds[fi++] = v3v0;	inds[fi++] = p[4];
				}
			}
		}break;

		case 4:
		{
		//	we only treat the two most important cases here: either all new
		//	vertices lie in the base or none of them does.
			if(cornerStatus[4] == 0){
			//	all lie in the base
			//	we require an additional face vertex here
				int nVrt = FACE_FROM_VRTS[0][1][2] + F;
				int v0v1 = EDGE_FROM_VRTS[0][1] + E;
				int v1v2 = EDGE_FROM_VRTS[1][2] + E;
				int v2v3 = EDGE_FROM_VRTS[2][3] + E;
				int v3v0 = EDGE_FROM_VRTS[3][0] + E;

			//	create 4 new pyramids
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = 5;
				inds[fi++] = 0;		inds[fi++] = v0v1;	inds[fi++] = nVrt;
				inds[fi++] = v3v0;	inds[fi++] = 4;

				inds[fi++] = 5;
				inds[fi++] = 1;		inds[fi++] = v1v2;	inds[fi++] = nVrt;
				inds[fi++] = v0v1;	inds[fi++] = 4;

				inds[fi++] = 5;
				inds[fi++] = 2;		inds[fi++] = v2v3;	inds[fi++] = nVrt;
				inds[fi++] = v1v2;	inds[fi++] = 4;

				inds[fi++] = 5;
				inds[fi++] = 3;		inds[fi++] = v3v0;	inds[fi++] = nVrt;
				inds[fi++] = v2v3;	inds[fi++] = 4;
			}
			else if(cornerStatus[4] == 4){
			//	all lie on the side edges
				int e0 = EDGE_FROM_VRTS[0][4] + E;
				int e1 = EDGE_FROM_VRTS[1][4] + E;
				int e2 = EDGE_FROM_VRTS[2][4] + E;
				int e3 = EDGE_FROM_VRTS[3][4] + E;

			//	create a hexahedron and a new pyramid
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = 8;
				inds[fi++] = 0;		inds[fi++] = 1;
				inds[fi++] = 2;		inds[fi++] = 3;
				inds[fi++] = e0;	inds[fi++] = e1;
				inds[fi++] = e2;	inds[fi++] = e3;

				inds[fi++] = 5;
				inds[fi++] = e0;	inds[fi++] = e1;	inds[fi++] = e2;
				inds[fi++] = e3;	inds[fi++] = 4;
			}
		}break;

		case 8:	// regular refine
		{
		//	we have to create 6 new pyramids and 4 tetrahedrons
			int& fi = fillCount;
			int* inds = newIndsOut;
			inds[fi++] = 5;
			inds[fi++] = 0;		inds[fi++] = E;
			inds[fi++] = F;		inds[fi++] = E + 3;		inds[fi++] = E + 4;

			inds[fi++] = 5;
			inds[fi++] = E;			inds[fi++] = 1;
			inds[fi++] = E + 1;		inds[fi++] = F;		inds[fi++] = E + 5;

			inds[fi++] = 5;
			inds[fi++] = F;		inds[fi++] = E + 1;
			inds[fi++] = 2;		inds[fi++] = E + 2;		inds[fi++] = E + 6;

			inds[fi++] = 5;
			inds[fi++] = E + 3;		inds[fi++] = F;
			inds[fi++] = E + 2;		inds[fi++] = 3;		inds[fi++] = E + 7;

			inds[fi++] = 5;
			inds[fi++] = E + 7;		inds[fi++] = E + 6;
			inds[fi++] = E + 5;		inds[fi++] = E + 4;		inds[fi++] = F;

			inds[fi++] = 5;
			inds[fi++] = E + 4;		inds[fi++] = E + 5;
			inds[fi++] = E + 6;		inds[fi++] = E + 7;		inds[fi++] = 4;

		//	tetrahedrons
			inds[fi++] = 4;
			inds[fi++] = E + 5;		inds[fi++] = E;
			inds[fi++] = E + 4;		inds[fi++] = F;

			inds[fi++] = 4;
			inds[fi++] = E + 6;		inds[fi++] = E + 1;
			inds[fi++] = E + 5;		inds[fi++] = F;

			inds[fi++] = 4;
			inds[fi++] = E + 7;		inds[fi++] = E + 2;
			inds[fi++] = E + 6;		inds[fi++] = F;

			inds[fi++] = 4;
			inds[fi++] = E + 4;		inds[fi++] = E + 3;
			inds[fi++] = E + 7;		inds[fi++] = F;

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

}// end of namespace
}// end of namespace

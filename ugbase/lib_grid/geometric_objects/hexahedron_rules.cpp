// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.05.2011 (m,d,y)

#include <cassert>
#include "rule_util.h"
#include "hexahedron_rules.h"

namespace ug{
namespace hex_rules{

///	Output are the vertices of a rotated quadrilateral.
void RotateQuad(int vrtsOut[4], const int quad[4], int steps)
{
	if(steps < 0)
		steps = (4 - ((-steps) % 4));

	for(int i = 0; i < 4; ++i)
		vrtsOut[(i + steps) % 4] = quad[i];
}

int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut, vector3*)
{
	newCenterOut = false;
//	If a refinement rule is not implemented, fillCount will stay at 0.
//	Before returning, we'll check, whether fillCount is at 0 and will
//	perform refinement with an additional vertex in this case.

//	corner status is used to mark vertices, which are connected to refined edges
//	the value tells, how many associated edges are refined.
	int cornerStatus[8] = {0, 0, 0, 0, 0, 0, 0, 0};

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
	const int V = NUM_VERTICES + NUM_EDGES + NUM_FACES;

//	depending on the number of new vertices, we will now apply different
//	refinement rules. Further distinction may have to be done.
	switch(numNewVrts){
		case 0:
		{
		//	simply put the default prism back to newIndsOut
			newIndsOut[fillCount++] = 8;
			newIndsOut[fillCount++] = 0;
			newIndsOut[fillCount++] = 1;
			newIndsOut[fillCount++] = 2;
			newIndsOut[fillCount++] = 3;
			newIndsOut[fillCount++] = 4;
			newIndsOut[fillCount++] = 5;
			newIndsOut[fillCount++] = 6;
			newIndsOut[fillCount++] = 7;
		}break;

		case 1:
		{
		//	we'll iterate over all faces. If a face does not contain the
		//	refined edge, then we'll create a pyramid from it.
			const int refEdge = refEdgeInds[0];
			const int nVrt = refEdge + E;

			int& fi = fillCount;
			int* inds = newIndsOut;

			for(int i = 0; i < NUM_FACES; ++i){
				if(!FACE_CONTAINS_EDGE[i][refEdge]){
					const int* f = FACE_VRT_INDS[i];
					inds[fi++] = 5;
					inds[fi++] = f[0];	inds[fi++] = f[1];
					inds[fi++] = f[2];	inds[fi++] = f[3];
					inds[fi++] = nVrt;
				}
			}
		}break;

		case 2:
		{
		//	currently only the case where both new vertices lie on opposed
		//	edges of one side is supported.
		//	We thus try to get that face.
			int refFaceInd = FACE_FROM_EDGES[refEdgeInds[0]][refEdgeInds[1]];
			if(refFaceInd == -1)
				break;

		//	make sure that corners of refEdge0 both have corner-state 1
			const int* refEdge0 = EDGE_VRT_INDS[refEdgeInds[0]];
			if(cornerStatus[refEdge0[0]] != 1 || cornerStatus[refEdge0[1]] != 1)
				break;

		//	find the first face != refFaceInd, which contains refEdge0
			int faceInd = -1;
			for(int i = 0; i < NUM_FACES; ++i){
				if(FACE_CONTAINS_EDGE[i][refEdgeInds[0]] && i != refFaceInd){
					faceInd = i;
					break;
				}
			}

			assert(faceInd != -1);

		//	the face and its opposing face
			const int* tf = FACE_VRT_INDS[faceInd];
			const int* tof = OPPOSED_FACE_VRT_INDS[faceInd];

		//	we have to rotate the faces, so that the rules can be easily applied
			int localEdgeInd = -1;
			for(int i = 0; i < 4; ++i){
				if(newEdgeVrts[EDGE_FROM_VRTS[tf[i]][tf[(i + 1) % 4]]]){
					localEdgeInd = i;
					break;
				}
			}

			assert(localEdgeInd != -1);

		//	the local refEdge shall lie between vrts 1, 2
			int f[4];
			int of[4];
			RotateQuad(f, tf, 1 - localEdgeInd);
			RotateQuad(of, tof, 1 - localEdgeInd);

		//	create three new prisms
			const int v0 = EDGE_FROM_VRTS[f[1]][f[2]] + E;
			const int v1 = EDGE_FROM_VRTS[of[1]][of[2]] + E;

			int& fi = fillCount;
			int* inds = newIndsOut;

			inds[fi++] = 6;
			inds[fi++] = f[0];	inds[fi++] = f[1];	inds[fi++] = v0;
			inds[fi++] = of[0];	inds[fi++] = of[1];	inds[fi++] = v1;

			inds[fi++] = 6;
			inds[fi++] = f[0];	inds[fi++] = v0;	inds[fi++] = f[3];
			inds[fi++] = of[0];	inds[fi++] = v1;	inds[fi++] = of[3];

			inds[fi++] = 6;
			inds[fi++] = f[2];	inds[fi++] = f[3];	inds[fi++] = v0;
			inds[fi++] = of[2];	inds[fi++] = of[3];	inds[fi++] = v1;
		}break;

		case 4:
		{
		//	we only support nice regular cuts with two resulting hexahedrons.
		//	get the face which contains refEdge0
			int faceInd = -1;
			for(int i = 0; i < NUM_FACES; ++i){
				if(FACE_CONTAINS_EDGE[i][refEdgeInds[0]]){
					faceInd = i;
					break;
				}
			}
			assert(faceInd != -1);

			const int* tf = FACE_VRT_INDS[faceInd];
			const int* tof = OPPOSED_FACE_VRT_INDS[faceInd];

		//	get the local index of refEdgeInds[0]
			int localEdgeInd = -1;
			for(int i = 0; i < 4; ++i){
				if(newEdgeVrts[EDGE_FROM_VRTS[tf[i]][tf[(i + 1) % 4]]]){
					localEdgeInd = i;
					break;
				}
			}
			assert(localEdgeInd != -1);

		//	Check whether the opposed edge in tf will be refined, too
			if(newEdgeVrts[EDGE_FROM_VRTS[tf[(localEdgeInd + 2) % 4]]
			                             [tf[(localEdgeInd + 3) % 4]]])
			{
			//	make sure that those are the same local indices as in tof
				if(!newEdgeVrts[EDGE_FROM_VRTS[tof[(localEdgeInd)]]
				                              [tof[(localEdgeInd + 1) % 4]]]
				   || !newEdgeVrts[EDGE_FROM_VRTS[tof[(localEdgeInd + 2) % 4]]
				                              [tof[(localEdgeInd + 3) % 4]]])
				{
				//	they are not. we can't handle that case.
					break;
				}
			}
			else
				break;

		//	we rotate the faces, so that the refined edges always lie
		//	between vertices 1,2 and 3,0
			int f[4];
			int of[4];

			RotateQuad(f, tf, 1 - localEdgeInd);
			RotateQuad(of, tof, 1 - localEdgeInd);

			int v0, v1, ov0, ov1;

		//	get the new vertex indices of new vertices
			v0 = EDGE_FROM_VRTS[f[1]][f[2]] + E;
			v1 = EDGE_FROM_VRTS[f[3]][f[0]] + E;
			ov0 = EDGE_FROM_VRTS[of[1]][of[2]] + E;
			ov1 = EDGE_FROM_VRTS[of[3]][of[0]] + E;

		//	create 2 new hexahedrons
			int& fi = fillCount;
			int* inds = newIndsOut;

			inds[fi++] = 8;
			inds[fi++] = f[0];	inds[fi++] = f[1];
			inds[fi++] = v0;	inds[fi++] = v1;
			inds[fi++] = of[0];	inds[fi++] = of[1];
			inds[fi++] = ov0;	inds[fi++] = ov1;

			inds[fi++] = 8;
			inds[fi++] = f[2];	inds[fi++] = f[3];
			inds[fi++] = v1;	inds[fi++] = v0;
			inds[fi++] = of[2];	inds[fi++] = of[3];
			inds[fi++] = ov1;	inds[fi++] = ov0;
		}break;

		case 8:
		{
		//	we currently only support the case where two opposing faces
		//	are completly refined.
		//	iterate through all faces, until two are found whose edges are all
		//	to be refined.
			int refFaceInd[2] = {-1, -1};
			for(int i = 0; i < NUM_FACES; ++i){
				const int* e = FACE_EDGE_INDS[i];
				int j = 0;
				for(;j < 4; ++j){
					if(!newEdgeVrts[e[j]])
						break;
				}

				if(j == 4){
					if(refFaceInd[0] == -1)
						refFaceInd[0] = i;
					else{
						assert(refFaceInd[1] == -1);
						refFaceInd[1] = i;
					}
				}
			}

		//	if two faces were found and if they are opposed faces, then
		//	we can create the four new hexahedrons
			if(refFaceInd[0] != -1 && refFaceInd[1] != -1
			   && refFaceInd[1] == OPPOSED_FACE[refFaceInd[0]])
			{
				const int* f = FACE_VRT_INDS[refFaceInd[0]];
				const int* of = OPPOSED_FACE_VRT_INDS[refFaceInd[0]];

				const int e0 = EDGE_FROM_VRTS[f[0]][f[1]] + E;
				const int e1 = EDGE_FROM_VRTS[f[1]][f[2]] + E;
				const int e2 = EDGE_FROM_VRTS[f[2]][f[3]] + E;
				const int e3 = EDGE_FROM_VRTS[f[3]][f[0]] + E;

				const int oe0 = EDGE_FROM_VRTS[of[0]][of[1]] + E;
				const int oe1 = EDGE_FROM_VRTS[of[1]][of[2]] + E;
				const int oe2 = EDGE_FROM_VRTS[of[2]][of[3]] + E;
				const int oe3 = EDGE_FROM_VRTS[of[3]][of[0]] + E;

				const int fvrt = refFaceInd[0] + F;
				const int ofvrt = refFaceInd[1] + F;

			//	create 4 new hexahedrons
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = 8;
				inds[fi++] = f[0];	inds[fi++] = e0;
				inds[fi++] = fvrt;	inds[fi++] = e3;
				inds[fi++] = of[0];	inds[fi++] = oe0;
				inds[fi++] = ofvrt;	inds[fi++] = oe3;

				inds[fi++] = 8;
				inds[fi++] = f[1];	inds[fi++] = e1;
				inds[fi++] = fvrt;	inds[fi++] = e0;
				inds[fi++] = of[1];	inds[fi++] = oe1;
				inds[fi++] = ofvrt;	inds[fi++] = oe0;

				inds[fi++] = 8;
				inds[fi++] = f[2];	inds[fi++] = e2;
				inds[fi++] = fvrt;	inds[fi++] = e1;
				inds[fi++] = of[2];	inds[fi++] = oe2;
				inds[fi++] = ofvrt;	inds[fi++] = oe1;

				inds[fi++] = 8;
				inds[fi++] = f[3];	inds[fi++] = e3;
				inds[fi++] = fvrt;	inds[fi++] = e2;
				inds[fi++] = of[3];	inds[fi++] = oe3;
				inds[fi++] = ofvrt;	inds[fi++] = oe2;
			}
		}break;

		case 12: // regular refine
		{
		//	we have to create 8 new hexahedrons
			int& fi = fillCount;
			int* inds = newIndsOut;
			inds[fi++] = 8;
			inds[fi++] = 0;			inds[fi++] = E;
			inds[fi++] = F;			inds[fi++] = E + 3;
			inds[fi++] = E + 4;		inds[fi++] = F + 1;
			inds[fi++] = V;			inds[fi++] = F + 4;

			inds[fi++] = 8;
			inds[fi++] = E;			inds[fi++] = 1;
			inds[fi++] = E + 1;		inds[fi++] = F;
			inds[fi++] = F + 1;		inds[fi++] = E + 5;
			inds[fi++] = F + 2;		inds[fi++] = V;

			inds[fi++] = 8;
			inds[fi++] = F;			inds[fi++] = E + 1;
			inds[fi++] = 2;			inds[fi++] = E + 2;
			inds[fi++] = V;			inds[fi++] = F + 2;
			inds[fi++] = E + 6;		inds[fi++] = F + 3;

			inds[fi++] = 8;
			inds[fi++] = E + 3;		inds[fi++] = F;
			inds[fi++] = E + 2;		inds[fi++] = 3;
			inds[fi++] = F + 4;		inds[fi++] = V;
			inds[fi++] = F + 3;		inds[fi++] = E + 7;

			inds[fi++] = 8;
			inds[fi++] = E + 4;		inds[fi++] = F + 1;
			inds[fi++] = V;			inds[fi++] = F + 4;
			inds[fi++] = 4;			inds[fi++] = E + 8;
			inds[fi++] = F + 5;		inds[fi++] = E + 11;

			inds[fi++] = 8;
			inds[fi++] = F + 1;		inds[fi++] = E + 5;
			inds[fi++] = F + 2;		inds[fi++] = V;
			inds[fi++] = E + 8;		inds[fi++] = 5;
			inds[fi++] = E + 9;		inds[fi++] = F + 5;

			inds[fi++] = 8;
			inds[fi++] = V;			inds[fi++] = F + 2;
			inds[fi++] = E + 6;		inds[fi++] = F + 3;
			inds[fi++] = F + 5;		inds[fi++] = E + 9;
			inds[fi++] = 6;			inds[fi++] = E + 10;

			inds[fi++] = 8;
			inds[fi++] = F + 4;		inds[fi++] = V;
			inds[fi++] = F + 3;		inds[fi++] = E + 7;
			inds[fi++] = E + 11;	inds[fi++] = F + 5;
			inds[fi++] = E + 10;	inds[fi++] = 7;

		//	the rule requires a new center vertex
			newCenterOut = true;
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

// created by mstepnie
// adapted from Sebastian Reiter
// martin.stepniewski@gcsc.uni-frankfurt.de
// 10.07.2014 (m,d,y)

#include <cassert>
#include "octahedron_rules.h"
#include "rule_util.h"
#include "grid_object_ids.h"

namespace ug{
namespace oct_rules
{

///	Output are the vertices of a octahedron rotated around its vertical axis
void RotateOctahedron(int vrtsOut[NUM_VERTICES], int steps)
{
	if(steps < 0)
		steps = (4 - ((-steps) % 4));

	for(int i = 0; i < 4; ++i)
		vrtsOut[((i + steps) % 4) + 1] = i+1;

	vrtsOut[0] = 0;
	vrtsOut[5] = 5;
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
	//const int F = NUM_VERTICES + NUM_EDGES;
	const int V = NUM_VERTICES + NUM_EDGES + NUM_FACES;

//	depending on the number of new vertices, we will now apply different
//	refinement rules. Further distinction may have to be done.
	switch(numNewVrts){
		case 0:
		{
		//	simply put the default octahedron back to newIndsOut
			newIndsOut[fillCount++] = GOID_OCTAHEDRON;
			newIndsOut[fillCount++] = 0;
			newIndsOut[fillCount++] = 1;
			newIndsOut[fillCount++] = 2;
			newIndsOut[fillCount++] = 3;
			newIndsOut[fillCount++] = 4;
			newIndsOut[fillCount++] = 5;
		}break;

	//	TODO: Adaptive refinement cases still to be implemented!

	//	REGULAR REFINEMENT
		case 12:
		{
		//	We have to create 6 octahedrons and 8 tetrahedrons
			int& fi = fillCount;
			int* inds = newIndsOut;

		//	octahedrons
			inds[fi++] = GOID_OCTAHEDRON;
			inds[fi++] = 0;			inds[fi++] = E;
			inds[fi++] = E + 1;		inds[fi++] = E + 2;
			inds[fi++] = E + 3;		inds[fi++] = V;

			inds[fi++] = GOID_OCTAHEDRON;
			inds[fi++] = E;			inds[fi++] = 1;
			inds[fi++] = E + 4;		inds[fi++] = V;
			inds[fi++] = E + 7;		inds[fi++] = E + 8;

			inds[fi++] = GOID_OCTAHEDRON;
			inds[fi++] = E + 1;		inds[fi++] = E + 4;
			inds[fi++] = 2;			inds[fi++] = E + 5;
			inds[fi++] = V;			inds[fi++] = E + 9;

			inds[fi++] = GOID_OCTAHEDRON;
			inds[fi++] = E + 2;		inds[fi++] = V;
			inds[fi++] = E + 5;		inds[fi++] = 3;
			inds[fi++] = E + 6;		inds[fi++] = E + 10;

			inds[fi++] = GOID_OCTAHEDRON;
			inds[fi++] = E + 3;		inds[fi++] = E + 7;
			inds[fi++] = V;			inds[fi++] = E + 6;
			inds[fi++] = 4;			inds[fi++] = E + 11;

			inds[fi++] = GOID_OCTAHEDRON;
			inds[fi++] = V;			inds[fi++] = E + 8;
			inds[fi++] = E + 9;		inds[fi++] = E + 10;
			inds[fi++] = E + 11;	inds[fi++] = 5;

		//	lower tetrahedrons
			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 4;		inds[fi++] = E + 1;
			inds[fi++] = E;			inds[fi++] = V;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 5;		inds[fi++] = E + 2;
			inds[fi++] = E + 1;		inds[fi++] = V;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 6;		inds[fi++] = E + 3;
			inds[fi++] = E + 2;		inds[fi++] = V;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 7;		inds[fi++] = E;
			inds[fi++] = E + 3;		inds[fi++] = V;

		//	upper tetrahedrons
			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 9;		inds[fi++] = E + 4;
			inds[fi++] = E + 8;		inds[fi++] = V;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 10;	inds[fi++] = E + 5;
			inds[fi++] = E + 9;		inds[fi++] = V;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 11;	inds[fi++] = E + 6;
			inds[fi++] = E + 10;	inds[fi++] = V;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = E + 8;		inds[fi++] = E + 7;
			inds[fi++] = E + 11;	inds[fi++] = V;

		//	the rule requires a new center vertex
			newCenterOut = true;

			/*
			 * Former orientation
			 *
		//	lower tetrahedrons
			inds[fi++] = 4;
			inds[fi++] = E;			inds[fi++] = E + 1;
			inds[fi++] = V;			inds[fi++] = E + 4;

			inds[fi++] = 4;
			inds[fi++] = E + 1;		inds[fi++] = E + 2;
			inds[fi++] = V;			inds[fi++] = E + 5;

			inds[fi++] = 4;
			inds[fi++] = E + 2;		inds[fi++] = E + 3;
			inds[fi++] = V;			inds[fi++] = E + 6;

			inds[fi++] = 4;
			inds[fi++] = E + 3;		inds[fi++] = E;
			inds[fi++] = V;			inds[fi++] = E + 7;

		//	upper tetrahedrons
			inds[fi++] = 4;
			inds[fi++] = E + 8;		inds[fi++] = E + 9;
			inds[fi++] = V;			inds[fi++] = E + 4;

			inds[fi++] = 4;
			inds[fi++] = E + 9;		inds[fi++] = E + 10;
			inds[fi++] = V;			inds[fi++] = E + 5;

			inds[fi++] = 4;
			inds[fi++] = E + 10;	inds[fi++] = E + 11;
			inds[fi++] = V;			inds[fi++] = E + 6;

			inds[fi++] = 4;
			inds[fi++] = E + 11;	inds[fi++] = E + 8;
			inds[fi++] = V;			inds[fi++] = E + 7;
			*/

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


}//	end of namespace tet_rules
}//	end of namespace ug

// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 31.05.2011 (m,d,y)
 
#include "hexahedron_rules.h"

namespace ug{
namespace hex_rules{


int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut)
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

	return fillCount;
}

}// end of namespace
}// end of namespace

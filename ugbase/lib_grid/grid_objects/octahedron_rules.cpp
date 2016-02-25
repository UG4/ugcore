/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Stepniewski
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


int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut, vector3* corners)
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
	// int refEdgeInds[NUM_EDGES];

//	count the number of new vertices and fill newVrtEdgeInds
	int numNewVrts = 0;
	for(int i = 0; i < NUM_EDGES; ++i){
		if(newEdgeVrts[i]){
			// refEdgeInds[numNewVrts] = i;
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

		//	For the 6 octahedrons we'll choose the shortest diagonals
		//	and order the octahedral nodes, so that, the implicit
		//	subdivision of the octahedron into tetrahedrons during
		//	discretization happens alongside this shortest diagonal.

			int bestDiag[6];

			for(int i = 0; i < 6; ++i)
				bestDiag[i] = 2;

			if(corners){
			//
			//	Calculate edge centers and element center
			//

				vector3 c06, c07, c08, c09;
				vector3 c10, c11, c12, c13;
				vector3 c14, c15, c16, c17;
				vector3 c18;

			//	Bottom hemisphere edges
				VecAdd(c06, corners[0], corners[1]);
				VecScale(c06, c06, 0.5);

				VecAdd(c07, corners[0], corners[2]);
				VecScale(c07, c07, 0.5);

				VecAdd(c08, corners[0], corners[3]);
				VecScale(c08, c08, 0.5);

				VecAdd(c09, corners[0], corners[4]);
				VecScale(c09, c09, 0.5);

			//	Middle ring edges
				VecAdd(c10, corners[1], corners[2]);
				VecScale(c10, c10, 0.5);

				VecAdd(c11, corners[2], corners[3]);
				VecScale(c11, c11, 0.5);

				VecAdd(c12, corners[3], corners[4]);
				VecScale(c12, c12, 0.5);

				VecAdd(c13, corners[4], corners[1]);
				VecScale(c13, c13, 0.5);

			//	Top hemisphere edges
				VecAdd(c14, corners[1], corners[5]);
				VecScale(c14, c14, 0.5);

				VecAdd(c15, corners[2], corners[5]);
				VecScale(c15, c15, 0.5);

				VecAdd(c16, corners[3], corners[5]);
				VecScale(c16, c16, 0.5);

				VecAdd(c17, corners[4], corners[5]);
				VecScale(c17, c17, 0.5);

			//	Element center
				VecAdd(c18, corners[0], corners[5]);
				VecAppend(c18, corners[1], corners[2], corners[3], corners[4]);
				VecScale(c18, c18, 1.0/6.0);


			//
			//	Calculate the three diagonals of all 6 sub-octahedrons
			//

			//	Sub-octahedron 0
				number oct0_d0_0608 = VecDistanceSq(c06, c08);
				number oct0_d1_0709 = VecDistanceSq(c07, c09);
				number oct0_d2_0018 = VecDistanceSq(corners[0], c18);

				number d = oct0_d2_0018;

				if(oct0_d1_0709 < d)
				{
					bestDiag[0] = 1;
					d = oct0_d1_0709;
				}

				if(oct0_d0_0608 < d)
				{
					bestDiag[0] = 0;
				}

				//UG_LOG("oct0_d0_0608: " << oct0_d0_0608 << "; \t oct0_d1_0709: " << oct0_d1_0709 << "; \t oct0_d2_0018: " << oct0_d2_0018 << "; \t bestDiag0: " << bestDiag[0] << std::endl);

			//	Sub-octahedron 1
				number oct1_d0_0118 = VecDistanceSq(corners[1], c18);
				number oct1_d1_1013 = VecDistanceSq(c10, c13);
				number oct1_d2_0614 = VecDistanceSq(c06, c14);

				d = oct1_d2_0614;

				if(oct1_d1_1013 < d)
				{
					bestDiag[1] = 1;
					d = oct1_d1_1013;
				}

				if(oct1_d0_0118 < d)
				{
					bestDiag[1] = 0;
				}

				//UG_LOG("oct1_d0_0118: " << oct1_d0_0118 << "; \t oct1_d1_1013: " << oct1_d1_1013 << ";  \t oct1_d2_0614: " << oct1_d2_0614 << "; \t bestDiag1: " << bestDiag[1] << std::endl);

			//	Sub-octahedron 2
				number oct2_d0_1011 = VecDistanceSq(c10, c11);
				number oct2_d1_0218 = VecDistanceSq(corners[2], c18);
				number oct2_d2_0715 = VecDistanceSq(c07, c15);

				d = oct2_d2_0715;

				if(oct2_d1_0218 < d)
				{
					bestDiag[2] = 1;
					d = oct2_d1_0218;
				}

				if(oct2_d0_1011 < d)
				{
					bestDiag[2] = 0;
				}

				//UG_LOG("oct2_d0_1011: " << oct2_d0_1011 << "; \t oct2_d1_0218: " << oct2_d1_0218 << "; \t oct2_d2_0715: " << oct2_d2_0715 << "; \t bestDiag2: " << bestDiag[2] << std::endl);

			//	Sub-octahedron 3
				number oct3_d0_1803 = VecDistanceSq(c18, corners[3]);
				number oct3_d1_1112 = VecDistanceSq(c11, c12);
				number oct3_d2_0816 = VecDistanceSq(c08, c16);

				d = oct3_d2_0816;

				if(oct3_d1_1112 < d)
				{
					bestDiag[3] = 1;
					d = oct3_d1_1112;
				}

				if(oct3_d0_1803 < d)
				{
					bestDiag[3] = 0;
				}

				//UG_LOG("oct3_d0_1803: " << oct3_d0_1803 << "; \t oct3_d1_1112: " << oct3_d1_1112 << "; \t oct3_d2_0816: " << oct3_d2_0816 << "; \t bestDiag3: " << bestDiag[3] << std::endl);

			//	Sub-octahedron 4
				number oct4_d0_1312 = VecDistanceSq(c13, c12);
				number oct4_d1_1804 = VecDistanceSq(c18, corners[4]);
				number oct4_d2_0917 = VecDistanceSq(c09, c17);

				d = oct4_d2_0917;

				if(oct4_d1_1804 < d)
				{
					bestDiag[4] = 1;
					d = oct4_d1_1804;
				}

				if(oct4_d0_1312 < d)
				{
					bestDiag[4] = 0;
				}

				//UG_LOG("oct4_d0_1312: " << oct4_d0_1312 << "; \t oct4_d1_1804: " << oct4_d1_1804 << "; \t oct4_d2_0917: " << oct4_d2_0917 << "; \t bestDiag4: " << bestDiag[4] << std::endl);


			//	Sub-octahedron 5
				number oct5_d0_1416 = VecDistanceSq(c14, c16);
				number oct5_d1_1517 = VecDistanceSq(c15, c17);
				number oct5_d2_1805 = VecDistanceSq(c18, corners[5]);

				d = oct5_d0_1416;

				if(oct5_d1_1517 < d)
				{
					bestDiag[5] = 1;
					d = oct5_d1_1517;
				}

				if(oct5_d2_1805 < d)
				{
					bestDiag[5] = 0;
				}

				//UG_LOG("oct5_d0_1416: " << oct5_d0_1416 << "; \t oct5_d1_1517: " << oct5_d1_1517 << "; \t oct5_d2_1805: " << oct5_d2_1805 << "; \t bestDiag5: " << bestDiag[5] << std::endl);
			}

		//	Experimental fixed choice of the bestDiag
			for(int i = 0; i < 6; ++i)
				bestDiag[i] = 0;

		//
		//	define octahedrons
		//

//			for(int i = 0; i < 6; ++i)
//				UG_LOG(bestDiag[i] << "; ");
//			UG_LOG(std::endl);

		//	Sub-octahedron 0
			switch(bestDiag[0]){
				case 0: // diag 6-8
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = 0;			inds[fi++] = E;
					inds[fi++] = E + 1;		inds[fi++] = E + 2;
					inds[fi++] = E + 3;		inds[fi++] = V;
				break;

				case 1: // diag 7-9
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = 0;			inds[fi++] = E + 3;
					inds[fi++] = E;			inds[fi++] = E + 1;
					inds[fi++] = E + 2;		inds[fi++] = V;
				break;

				case 2: // diag 0-18
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E;			inds[fi++] = V;
					inds[fi++] = E + 1;		inds[fi++] = 0;
					inds[fi++] = E + 3;		inds[fi++] = E + 2;
				break;
			}

		//	Sub-octahedron 1
			switch(bestDiag[1]){
				case 0: // diag 1-18
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E;			inds[fi++] = 1;
					inds[fi++] = E + 4;		inds[fi++] = V;
					inds[fi++] = E + 7;		inds[fi++] = E + 8;
				break;

				case 1: // diag 10-13
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E;			inds[fi++] = E + 7;
					inds[fi++] = 1;			inds[fi++] = E + 4;
					inds[fi++] = V;			inds[fi++] = E + 8;
				break;

				case 2: // diag 6-14
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = 1;			inds[fi++] = E + 8;
					inds[fi++] = E + 4;		inds[fi++] = E;
					inds[fi++] = E + 7;		inds[fi++] = V;
				break;
			}

		//	Sub-octahedron 2
			switch(bestDiag[2]){
				case 0: // diag 10-11
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 1;		inds[fi++] = E + 4;
					inds[fi++] = 2;			inds[fi++] = E + 5;
					inds[fi++] = V;			inds[fi++] = E + 9;
				break;

				case 1: // diag 2-18
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 1;		inds[fi++] = V;
					inds[fi++] = E + 4;		inds[fi++] = 2;
					inds[fi++] = E + 5;		inds[fi++] = E + 9;
				break;

				case 2: // diag 7-15
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 4;		inds[fi++] = E + 9;
					inds[fi++] = 2; 		inds[fi++] = E + 1;
					inds[fi++] = V; 		inds[fi++] = E + 5;
				break;
			}

		//	Sub-octahedron 3
			switch(bestDiag[3]){
				case 0: // diag 18-3
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 2;		inds[fi++] = V;
					inds[fi++] = E + 5;		inds[fi++] = 3;
					inds[fi++] = E + 6;		inds[fi++] = E + 10;
				break;

				case 1: // diag 11-12
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 2;		inds[fi++] = E + 6;
					inds[fi++] = V;		    inds[fi++] = E + 5;
					inds[fi++] = 3;			inds[fi++] = E + 10;
				break;

				case 2: // diag 8-16
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = V;			inds[fi++] = E + 10;
					inds[fi++] = E + 5;		inds[fi++] = E + 2;
					inds[fi++] = E + 6;		inds[fi++] = 3;
				break;
			}

		//	Sub-octahedron 4
			switch(bestDiag[4]){
				case 0: // diag 13-12
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 3;		inds[fi++] = E + 7;
					inds[fi++] = V;			inds[fi++] = E + 6;
					inds[fi++] = 4;			inds[fi++] = E + 11;
				break;

				case 1: // diag 18-4
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 3;		inds[fi++] = 4;
					inds[fi++] = E + 7;		inds[fi++] = V;
					inds[fi++] = E + 6;		inds[fi++] = E + 11;
				break;

				case 2: // diag 9-17
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 7;		inds[fi++] = E + 11;
					inds[fi++] = V;			inds[fi++] = E + 3;
					inds[fi++] = 4;			inds[fi++] = E + 6;
				break;
			}

		//	Sub-octahedron 5
			switch(bestDiag[5]){
				case 0: // diag 14-16
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = V;			inds[fi++] = E + 8;
					inds[fi++] = E + 9;		inds[fi++] = E + 10;
					inds[fi++] = E + 11;	inds[fi++] = 5;
				break;

				case 1: // diag 15-17
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = V;			inds[fi++] = E + 11;
					inds[fi++] = E + 8;		inds[fi++] = E + 9;
					inds[fi++] = E + 10;	inds[fi++] = 5;
				break;

				case 2: // diag 18-5
					inds[fi++] = GOID_OCTAHEDRON;
					inds[fi++] = E + 8;		inds[fi++] = 5;
					inds[fi++] = E + 9;		inds[fi++] = V;
					inds[fi++] = E + 11;	inds[fi++] = E + 10;
				break;
			}

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


}//	end of namespace oct_rules
}//	end of namespace ug

#include <cassert>
#include "tetrahedron_rules.h"
#include "rule_util.h"
#include "grid_object_ids.h"

namespace ug{
namespace tet_rules
{

/// global refinement rule information switching between regular and subdivision volume refinement
static GlobalRefinementRule g_refinementRule = STANDARD;

void SetRefinementRule(GlobalRefinementRule refRule)
{
	g_refinementRule = refRule;
}

GlobalRefinementRule GetRefinementRule()
{
	return g_refinementRule;
}


///	Rotates the given tetrahedron while keeping the specified point fixed
/**	The face opposed to the fixed point is rotated in counter-clockwise order when
 * viewed from the fixed point.
 * e.g.
 * \code
 *		int vrts[] = {0, 1, 2, 3};
 *		TetRotation(vrts, 3, 1);
 *	//	-> vrts == {2, 0, 1, 3}
 * \endcode
 */
void TetRotation (
		int vrtsInOut[NUM_VERTICES],
		const int fixedPoint,
		const size_t steps)
{
//	get the opposing face of the fixed point
	const int opFace = OPPOSED_OBJECT[fixedPoint][1];
	const int* finds = FACE_VRT_INDS[opFace];

	for(size_t i = 0; i < steps; ++i){
		const int tmp = vrtsInOut[finds[0]];
		vrtsInOut[finds[0]] = vrtsInOut[finds[2]];
		vrtsInOut[finds[2]] = vrtsInOut[finds[1]];
		vrtsInOut[finds[1]] = tmp;
	}
}

void InverseTetTransform(int* indsOut, const int* transformedInds){
	UG_ASSERT(indsOut != transformedInds, "The arrays have to differ!");
	for(int i = 0; i < NUM_VERTICES; ++i)
		indsOut[transformedInds[i]] = i;
}


int Refine(int* newIndsOut, int* newEdgeVrts, bool& newCenterOut,
		   vector3* corners)
{
	newCenterOut = false;
//	If a refinement rule is not implemented, fillCount will stay at 0.
//	Before returning, we'll check, whether fillCount is at 0 and will
//	perform refinement with an additional vertex in this case.

//	corner status is used to mark vertices, which are connected to refined edges
//	the value tells, how many associated edges are refined.
	int cornerStatus[4] = {0, 0, 0, 0};

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

//	depending on the number of new vertices, we will now apply different
//	refinement rules. Further distinction may have to be done.
	switch(numNewVrts){
		case 0:
		{
		//	simply put the default tetrahedron back to newIndsOut
			newIndsOut[fillCount++] = GOID_TETRAHEDRON;
			newIndsOut[fillCount++] = 0;
			newIndsOut[fillCount++] = 1;
			newIndsOut[fillCount++] = 2;
			newIndsOut[fillCount++] = 3;
		}break;

		case 1:
		{
			int refEdgeInd = refEdgeInds[0];
		//	the two faces which are not connected to the refined edge
		//	serve as bottom for the new tetrahedrons.
			for(int i_face = 0; i_face < NUM_FACES; ++i_face){
				if(!FACE_CONTAINS_EDGE[i_face][refEdgeInd]){
					const int* fvi = FACE_VRT_INDS[i_face];

					newIndsOut[fillCount++] = GOID_TETRAHEDRON;
					newIndsOut[fillCount++] = fvi[0];
					newIndsOut[fillCount++] = fvi[1];
					newIndsOut[fillCount++] = fvi[2];
					newIndsOut[fillCount++] = NUM_VERTICES + refEdgeInd;
				}
			}
		}break;

		case 2:
		{
		//	two types exist: The two refined edges share a vertex or not.
			if(OPPOSED_EDGE[refEdgeInds[0]] == refEdgeInds[1]){
			//	they do not share an edge
			//	we create a local order from refEdgeInds[0], which
			//	always has to be an edge in the base-triangle and
			//	refEdgeInds[1], which always connects a base-corner with
			//	the top.
				const int v0 = EDGE_VRT_INDS[refEdgeInds[0]][0];
				const int v1 = EDGE_VRT_INDS[refEdgeInds[0]][1];
				const int v2 = EDGE_VRT_INDS[refEdgeInds[1]][0];
				const int v3 = EDGE_VRT_INDS[refEdgeInds[1]][1];
				const int n0 = refEdgeInds[0] + NUM_VERTICES;
				const int n1 = refEdgeInds[1] + NUM_VERTICES;

			//	from this local order we can now construct the 4 new tetrahedrons.
				int& fi = fillCount;
				int* inds = newIndsOut;
				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v0;	inds[fi++] = n0;
				inds[fi++] = v2; 	inds[fi++] = n1;

				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v0;	inds[fi++] = n0;
				inds[fi++] = n1; 	inds[fi++] = v3;

				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = n0;	inds[fi++] = v1;
				inds[fi++] = v2; 	inds[fi++] = n1;

				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = n0;	inds[fi++] = v1;
				inds[fi++] = n1; 	inds[fi++] = v3;
			}
			else{
			//	they share an edge
			//	We have to create a pyramid and a tetrahedron.
			//	get the triangle, which contains both refEdges
				int refTriInd = FACE_FROM_EDGES[refEdgeInds[0]][refEdgeInds[1]];
				const int* f = FACE_VRT_INDS[refTriInd];

			//	find the edge (v0, v1) in refTri, which won't be refined
				int v0 = -1; int v1 = -1; int v2 = -1;
				for(int i = 0; i < 3; ++i){
					v0 = f[i];
					v1 = f[(i+1)%3];
					v2 = f[(i+2)%3];
					if(cornerStatus[v0] == 1 && cornerStatus[v1] == 1)
						break;
				}

			//	make sure that v2 is connected to two refined edges
				assert(cornerStatus[v2] == 2);

			//	get the new vertex on edge v1v2 and on v2v0
				int v1v2 = EDGE_FROM_VRTS[v1][v2] + NUM_VERTICES;
				int v2v0 = EDGE_FROM_VRTS[v2][v0] + NUM_VERTICES;

			//	get the top (vertex with cornerState 0)
				int vtop = -1;
				for(int i = 0; i < NUM_VERTICES; ++i){
					if(cornerStatus[i] == 0){
						vtop = i;
						break;
					}
				}

			//	now lets build the pyramid
				int& fi = fillCount;
				int* inds = newIndsOut;
				inds[fi++] = GOID_PYRAMID;
				inds[fi++] = v0;	inds[fi++] = v1;
				inds[fi++] = v1v2; 	inds[fi++] = v2v0;
				inds[fi++] = vtop;

			//	and now the terahedron
				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v2;	inds[fi++] = vtop;
				inds[fi++] = v2v0;	inds[fi++] = v1v2;
			}
		}break;

		case 3:
		{
		//	different possibilities exist. First we'll treat the one,
		//	where all new vertices lie on the edges of one triangle.
		//	Note that refTriInd could be -1 after the call.
			int refTriInd = FACE_FROM_EDGES[refEdgeInds[0]][refEdgeInds[1]];
			if(refTriInd == FACE_FROM_EDGES[refEdgeInds[1]][refEdgeInds[2]])
			{
			//	all three lie in one plane (refTriInd has to be != -1 to get here)
			//	get the top (vertex with cornerState 0)
				int vtop = -1;
				for(int i = 0; i < NUM_VERTICES; ++i){
					if(cornerStatus[i] == 0){
						vtop = i;
						break;
					}
				}

				const int* f = FACE_VRT_INDS[refTriInd];

			//	create four new tetrahedrons
				const int v0 = f[0]; const int v1 = f[1]; const int v2 = f[2];
				const int v0v1 = EDGE_FROM_VRTS[v0][v1] + NUM_VERTICES;
				const int v1v2 = EDGE_FROM_VRTS[v1][v2] + NUM_VERTICES;
				const int v2v0 = EDGE_FROM_VRTS[v2][v0] + NUM_VERTICES;

				int& fi = fillCount;
				int* inds = newIndsOut;
				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v0;	inds[fi++] = vtop;
				inds[fi++] = v0v1;	inds[fi++] = v2v0;
				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v1;	inds[fi++] = vtop;
				inds[fi++] = v1v2;	inds[fi++] = v0v1;
				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v2;	inds[fi++] = vtop;
				inds[fi++] = v2v0;	inds[fi++] = v1v2;
				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v0v1;	inds[fi++] = vtop;
				inds[fi++] = v1v2;	inds[fi++] = v2v0;
			}
			else{
			//	we have to further distinguish.
			//	gather the corners with corner-state 3 and corner state 1
				int corner3 = -1;
				int freeCorner[NUM_VERTICES];
				int freeCount = 0;
				for(int i = 0; i < NUM_VERTICES; ++i){
					if(cornerStatus[i] == 3)
						corner3 = i;
					else
						freeCorner[freeCount++] = i;
				}

				if(corner3 != -1){
				//	a corner with corner state 3 exists.
					assert(freeCount == 3);

				// get the face which won't be refined (required for correct orientation)
					int freeTriInd = FACE_FROM_VRTS[freeCorner[0]][freeCorner[1]]
												   [freeCorner[2]];

				//	the free tri is the new base, corner3 the top
					const int* f = FACE_VRT_INDS[freeTriInd];
					int v2v3 = EDGE_FROM_VRTS[f[2]][corner3] + NUM_VERTICES;
					int v1v3 = EDGE_FROM_VRTS[f[1]][corner3] + NUM_VERTICES;
					int v0v3 = EDGE_FROM_VRTS[f[0]][corner3] + NUM_VERTICES;

				//	add a prism and a tetrahedron
					int& fi = fillCount;
					int* inds = newIndsOut;
					inds[fi++] = GOID_PRISM;
					inds[fi++] = f[0]; inds[fi++] = f[1]; inds[fi++] = f[2];
					inds[fi++] = v0v3; inds[fi++] = v1v3; inds[fi++] = v2v3;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = v2v3; inds[fi++] = corner3; inds[fi++] = v0v3;
					inds[fi++] = v1v3;
				}
			}
		}break;

		case 4:
		{
		//	multiple settings with 4 refined edges exist.
		//	Either two opposing edges are not marked (case all2 == true)
		//	or the two unmarked edges are contained in 1 triangle

		//	check whether all vertices have corner state 2
			bool all2 = true;
			for(int i = 0; i < NUM_VERTICES; ++i){
				if(cornerStatus[i] !=2){
					all2 = false;
					break;
				}
			}

			if(all2){
			//	we've got a straight cut.
			//	we'll rotate the tetrahedron around the tip so, that edge 2 won't be refined.
				int steps = 0;
				if(!newEdgeVrts[EDGE_FROM_VRTS[0][1]])
					steps = 2;
				else if(!newEdgeVrts[EDGE_FROM_VRTS[1][2]])
					steps = 1;

				int t[NUM_VERTICES] = {0, 1, 2, 3};
				TetRotation(t, 3, steps);

				const int v0v1 = EDGE_FROM_VRTS[t[0]][t[1]] + NUM_VERTICES;
				const int v1v2 = EDGE_FROM_VRTS[t[1]][t[2]] + NUM_VERTICES;
				const int v0v3 = EDGE_FROM_VRTS[t[0]][t[3]] + NUM_VERTICES;
				const int v2v3 = EDGE_FROM_VRTS[t[2]][t[3]] + NUM_VERTICES;

				assert(newEdgeVrts[v0v1 - NUM_VERTICES]);
				assert(newEdgeVrts[v1v2 - NUM_VERTICES]);
				assert(newEdgeVrts[v0v3 - NUM_VERTICES]);
				assert(newEdgeVrts[v2v3 - NUM_VERTICES]);

			//	now build two prisms
				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = t[0];	inds[fi++] = v0v3;	inds[fi++] = v0v1;
				inds[fi++] = t[2];	inds[fi++] = v2v3;	inds[fi++] = v1v2;

				inds[fi++] = GOID_PRISM;
				inds[fi++] = t[1];	inds[fi++] = v1v2;	inds[fi++] = v0v1;
				inds[fi++] = t[3];	inds[fi++] = v2v3;	inds[fi++] = v0v3;
			}
			else{
			//	one corner has state 1, one has state 3 and two have state 2.
			//	Rotate the tet so that 1 is at the top and 4 is at the origin
				int I[NUM_VERTICES] = {0, 1, 2, 3};

				int s1 = -1;
				for(int i = 0; i < 4; ++i){
					if(cornerStatus[i] == 1){
						s1 = i;
						break;
					}
				}

				if(s1 != 3){
					const int fixedPoint = (s1 + 2) % 3;
					TetRotation(I, fixedPoint, 1);
				}

				if(cornerStatus[I[1]] == 3)
					TetRotation(I, 3, 2);
				else if(cornerStatus[I[2]] == 3)
					TetRotation(I, 3, 1);

			//	indices of edge-center vertices
				const int v01 = EDGE_FROM_VRTS[I[0]][I[1]] + NUM_VERTICES;
				const int v02 = EDGE_FROM_VRTS[I[0]][I[2]] + NUM_VERTICES;
				const int v03 = EDGE_FROM_VRTS[I[0]][I[3]] + NUM_VERTICES;
				const int v12 = EDGE_FROM_VRTS[I[1]][I[2]] + NUM_VERTICES;

				int& fi = fillCount;
				int* inds = newIndsOut;

				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = I[0];	inds[fi++] = v01;
				inds[fi++] = v02;	inds[fi++] = v03;

				inds[fi++] = GOID_TETRAHEDRON;
				inds[fi++] = v01;	inds[fi++] = v12;
				inds[fi++] = v02;	inds[fi++] = v03;

				inds[fi++] = GOID_PYRAMID;
				inds[fi++] = I[3];	inds[fi++] = I[1];
				inds[fi++] = v01;	inds[fi++] = v03;
				inds[fi++] = v12;

				inds[fi++] = GOID_PYRAMID;
				inds[fi++] = I[2];	inds[fi++] = I[3];
				inds[fi++] = v03;	inds[fi++] = v02;
				inds[fi++] = v12;
			}
		}break;

		case 5:{
		//	only one edge is not marked for refinement
			int unmarkedEdge = 0;
			for(int i = 0; i < NUM_EDGES; ++i){
				if(!newEdgeVrts[i]){
					unmarkedEdge = i;
					break;
				}
			}

			int uei0 = EDGE_VRT_INDS[unmarkedEdge][0];
			int uei1 = EDGE_VRT_INDS[unmarkedEdge][1];

		//	orientate the tetrahedron so that the unmarked edge connects
		//	vertices 2 and 3
			int I[NUM_VERTICES] = {0, 1, 2, 3};
		//	if the unmarked edge lies in the base triangle, we'll first have to
		//	rotate so that it connects a base-vertex with the top
			if(unmarkedEdge < 3){
				const int fixedPoint = (uei1 + 1) % 3;
				TetRotation(I, fixedPoint, 1);
				int IInv[4];
				InverseTetTransform(IInv, I);

				uei0 = IInv[uei0];
				uei1 = IInv[uei1];
				unmarkedEdge = EDGE_FROM_VRTS[uei0][uei1];
			}

		//	now rotate the tet to its final place
			if(unmarkedEdge == 3)
				TetRotation(I, 3, 2);
			else if(unmarkedEdge == 4)
				TetRotation(I, 3, 1);

		//	We obtained the final permutation I. Now create new elements
		//	indices of edge-center vertices
			const int v01 = EDGE_FROM_VRTS[I[0]][I[1]] + NUM_VERTICES;
			const int v02 = EDGE_FROM_VRTS[I[0]][I[2]] + NUM_VERTICES;
			const int v03 = EDGE_FROM_VRTS[I[0]][I[3]] + NUM_VERTICES;
			const int v12 = EDGE_FROM_VRTS[I[1]][I[2]] + NUM_VERTICES;
			const int v13 = EDGE_FROM_VRTS[I[1]][I[3]] + NUM_VERTICES;

			int& fi = fillCount;
			int* inds = newIndsOut;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = I[0];	inds[fi++] = v01;
			inds[fi++] = v02;	inds[fi++] = v03;

			inds[fi++] = GOID_TETRAHEDRON;
			inds[fi++] = I[1];	inds[fi++] = v12;
			inds[fi++] = v01;	inds[fi++] = v13;

			inds[fi++] = GOID_PYRAMID;
			inds[fi++] = v02;	inds[fi++] = v12;
			inds[fi++] = v13;	inds[fi++] = v03;
			inds[fi++] = v01;

			inds[fi++] = GOID_PRISM;
			inds[fi++] = v02;	inds[fi++] = v12;	inds[fi++] = I[2];
			inds[fi++] = v03;	inds[fi++] = v13;	inds[fi++] = I[3];


		}break;

	//	REGULAR REFINEMENT
		case 6:
		{
		//	depending on the user specified global refinement rule, we will now apply different
		//	refinement rules.
			switch(g_refinementRule){
				case STANDARD:
				{
					int& fi = fillCount;
					int* inds = newIndsOut;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = 0;					inds[fi++] = NUM_VERTICES + 0;
					inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 3;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = 1;					inds[fi++] = NUM_VERTICES + 1;
					inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 4;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = 2;					inds[fi++] = NUM_VERTICES + 2;
					inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 5;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 4;
					inds[fi++] = NUM_VERTICES + 5;	inds[fi++] = 3;

				//	for the remaining four tetrahedrons, we'll choose the shortest diagonal
					int bestDiag = 2;
					if(corners){
					//	there are three diagonals between the following edge-centers:
					//	0-5, 1-3, 2-4
						vector3 c1, c2;

					//	0-5
						VecAdd(c1, corners[0], corners[1]);
						VecScale(c1, c1, 0.5);
						VecAdd(c2, corners[2], corners[3]);
						VecScale(c2, c2, 0.5);
						number d05 = VecDistanceSq(c1, c2);

					//	1-3
						VecAdd(c1, corners[1], corners[2]);
						VecScale(c1, c1, 0.5);
						VecAdd(c2, corners[0], corners[3]);
						VecScale(c2, c2, 0.5);
						number d13 = VecDistanceSq(c1, c2);

					//	2-4
						VecAdd(c1, corners[0], corners[2]);
						VecScale(c1, c1, 0.5);
						VecAdd(c2, corners[1], corners[3]);
						VecScale(c2, c2, 0.5);
						number d = VecDistanceSq(c1, c2);

						if(d13 < d){
							bestDiag = 1;
							d = d13;
						}
						if(d05 < d){
							bestDiag = 0;
						}
					}

					switch(bestDiag){
					case 0:// diag: 0-5
						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 1;
						inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 4;
						inds[fi++] = NUM_VERTICES + 5;	inds[fi++] = NUM_VERTICES + 0;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;
						inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 0;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 3;
						inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 5;

						break;

					case 1:// diag: 1-3
						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 1;
						inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 3;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 4;
						inds[fi++] = NUM_VERTICES + 5;	inds[fi++] = NUM_VERTICES + 3;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;
						inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 1;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 3;
						inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 1;

						break;

					case 2:// diag 2-4
						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 2;
						inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 4;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 2;
						inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 4;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 3;
						inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 5;

						inds[fi++] = GOID_TETRAHEDRON;
						inds[fi++] = NUM_VERTICES + 4;	inds[fi++] = NUM_VERTICES + 1;
						inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 5;

						break;
					}
				}break;
				
				case HYBRID_TET_OCT:
				{
				//	REGULAR REFINEMENT
					int& fi = fillCount;
					int* inds = newIndsOut;

				//	outer tetrahedrons
					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = 0;					inds[fi++] = NUM_VERTICES + 0;
					inds[fi++] = NUM_VERTICES + 2;	inds[fi++] = NUM_VERTICES + 3;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = 1;					inds[fi++] = NUM_VERTICES + 1;
					inds[fi++] = NUM_VERTICES + 0;	inds[fi++] = NUM_VERTICES + 4;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = 2;					inds[fi++] = NUM_VERTICES + 2;
					inds[fi++] = NUM_VERTICES + 1;	inds[fi++] = NUM_VERTICES + 5;

					inds[fi++] = GOID_TETRAHEDRON;
					inds[fi++] = NUM_VERTICES + 3;	inds[fi++] = NUM_VERTICES + 4;
					inds[fi++] = NUM_VERTICES + 5;	inds[fi++] = 3;

				//	inner octahedron
					inds[fi++] = GOID_OCTAHEDRON;
					/*
					inds[fi++] = NUM_VERTICES + 1;		inds[fi++] = NUM_VERTICES + 0;
					inds[fi++] = NUM_VERTICES + 4;		inds[fi++] = NUM_VERTICES + 5;
					inds[fi++] = NUM_VERTICES + 2;		inds[fi++] = NUM_VERTICES + 3;
					*/
					inds[fi++] = NUM_VERTICES + 0;		inds[fi++] = NUM_VERTICES + 3;
					inds[fi++] = NUM_VERTICES + 4;		inds[fi++] = NUM_VERTICES + 1;
					inds[fi++] = NUM_VERTICES + 2;		inds[fi++] = NUM_VERTICES + 5;

				}break;
			}
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

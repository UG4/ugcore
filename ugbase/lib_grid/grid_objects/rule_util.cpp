// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 01.06.2011 (m,d,y)

#include <cassert>
#include <algorithm>
#include <iostream>
#include "rule_util.h"
#include "tetrahedron_rules.h"
#include "octahedron_rules.h"
#include "pyramid_rules.h"
#include "hexahedron_rules.h"

using namespace std;

namespace ug{
namespace shared_rules{

int RecursiveRefine(int* newIndsOut, int* newEdgeVrts,
					const int faceVrtInds[][4],
					const int faceEdgeInds[][4],
					int numVrts, int numEdges, int numFaces)
{
//	no refinement rule was found. We'll thus insert a new vertex in the
//	center and split the element into sub-elements.
//	we'll then recursively call this method for each sub-tetrahedron.
//	since all refinement rules for one side are given, one recursion will
//	always suffice.

//	a helper which avoids multiple recursion
	/*static bool dbgRecursionActive = false;
	assert(!dbgRecursionActive);
	dbgRecursionActive = true;*/


	int fillCount = 0;

//	we use constants that are chosen big enough to fit.
	int tmpInds[pyra_rules::MAX_NUM_INDS_OUT];
	int tmpNewEdgeVrts[12];

//	the indMap is used to map the generated temporary indices to the
//	indices of the main element.
	const int indMapSize = numVrts + numEdges + 1; // +1 for a quad-face
	int indMap[hex_rules::NUM_VERTICES + hex_rules::NUM_EDGES + 1];
	assert(indMapSize <= hex_rules::NUM_VERTICES + hex_rules::NUM_EDGES + 1);

	for(int i_face = 0; i_face < numFaces; ++i_face){
	//	fill indMap and tmpNewEdgeVrts.
		const int* f = faceVrtInds[i_face];
		int numVrtsOfNewElem = 5;
		if(f[3] == -1)
			numVrtsOfNewElem = 4;

		for(int i_ind = 0; i_ind < indMapSize; ++i_ind)
			indMap[i_ind] = -1;

		for(int i = 0; i < 4; ++i){
			if(f[i] != -1)
				indMap[i] = f[i];
		}

		for(int i_edge = 0; i_edge < numEdges; ++i_edge)
			tmpNewEdgeVrts[i_edge] = 0;

		for(int i_edge = 0; i_edge < 4; ++i_edge){
			const int edgeInd = faceEdgeInds[i_face][i_edge];
			if(edgeInd != -1 && newEdgeVrts[edgeInd]){
				tmpNewEdgeVrts[i_edge] = 1;
				indMap[numVrtsOfNewElem + i_edge] = edgeInd + numVrts;
			}
		}

	//	assign the index of the face vertex (only relevant for quad-faces) and
	//	assign index of inner vertex to indMap and call refine.
		bool centerVrtRequired = false;
		int count = 0;
		if(f[3] == -1){
			//clog << "calling tet_rules::Refine\n";
			indMap[3] = numVrts + numEdges + numFaces;
			count = tet_rules::Refine(tmpInds, tmpNewEdgeVrts, centerVrtRequired);
		}
		else{
			indMap[4] = numVrts + numEdges + numFaces;
			indMap[pyra_rules::NUM_VERTICES + pyra_rules::NUM_EDGES] =
													numVrts + numEdges + i_face;
			//clog << "calling pyra_rules::Refine\n";
			count = pyra_rules::Refine(tmpInds, tmpNewEdgeVrts, centerVrtRequired);
		}

	//	No center vertex should be created during this recudsion.
		assert(!centerVrtRequired);
/*
		clog << "tmpInds:";
		for(int i = 0; i < count; ++i)
			clog << " " << tmpInds[i];
		clog << endl;

		clog << "indMap:";
		for(int i = 0; i < indMapSize; ++i)
			clog << " " << indMap[i];
		clog << endl;
*/
	//	copy the new indices to the output array
		for(int i = 0; i < count;){
			int numElemInds = tmpInds[i];
			newIndsOut[fillCount + i] = numElemInds;
			++i;

			for(int j = 0; j < numElemInds; ++j, ++i){
				assert(indMap[tmpInds[i]] != -1);
				newIndsOut[fillCount + i] = indMap[tmpInds[i]];
			}
		}

		fillCount += count;
	}

	//dbgRecursionActive = false;
	return fillCount;
}

}// end of namespace
}// end of namespace

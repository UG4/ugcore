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
//	we'll then recursively call this method for each sub-tetrahedron/-pyramid.
//	since all refinement rules for one side exist for those element types,
//	one recursion will always suffice.

//	a helper which avoids multiple recursion
	/*static bool dbgRecursionActive = false;
	assert(!dbgRecursionActive);
	dbgRecursionActive = true;*/

	// #define LIBGRID_RECURSIVE_REFINE__PRINT_DEBUG_INFO
	#ifdef LIBGRID_RECURSIVE_REFINE__PRINT_DEBUG_INFO
		static bool firstVisit = true;
		if(firstVisit){
			firstVisit = false;
			UG_LOG("DEBUGGING: LIBGRID_RECURSIVE_REFINE__PRINT_DEBUG_INFO enabled.\n");
		}
		UG_LOG("DEBUG: RecursiveRefine called on element with " << numVrts << " vertices.\n");
		UG_LOG("       Edge mark pattern: ");
		for(int i = 0; i < numEdges; ++i){
			cout << " " << newEdgeVrts[i] != nullptr;
		}
		cout << endl;
	#endif


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

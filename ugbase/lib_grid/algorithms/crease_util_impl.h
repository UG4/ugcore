/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_crease_util_impl
#define __H__UG_crease_util_impl

#include <algorithm>
#include "lib_grid/iterators/associated_elements_iterator.h"

namespace ug{

template <class TVrtIter, class TAAPos>
void SelectKinkVertices(Selector& sel, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
						number thresholdAngle, bool selectDarts, TAAPos aaPos,
						Grid::edge_traits::callback cbConsiderEdge)
{
	using std::max;
	using std::endl;
	using vector_t = typename TAAPos::ValueType;

	if(!sel.grid())
		return;

	Grid& grid = *sel.grid();

	number thresholdRad = deg_to_rad(thresholdAngle);
	AssocElemIter<Vertex, Edge> eiter(cbConsiderEdge);

	for(TVrtIter viter = vrtsBegin; viter != vrtsEnd; ++viter){
		Vertex* v = *viter;
		Vertex* cv[2];
		int numCons = 0;
		for(eiter.reinit(grid, v); eiter.valid(); ++eiter){
			Edge* e = *eiter;
			if(numCons < 2)
				cv[numCons] = GetConnectedVertex(e, v);
			++numCons;
		}
		if(numCons == 2){
			vector_t dir0, dir1;
			VecSubtract(dir0, aaPos[v], aaPos[cv[0]]);
			VecSubtract(dir1, aaPos[cv[1]], aaPos[v]);
			if(VecAngle(dir0, dir1) >= thresholdRad)
				sel.select(v);
		}
		else if(selectDarts){
			sel.select(v);
		}
	}
}

}//	end of namespace

#endif	//__H__UG_crease_util_impl

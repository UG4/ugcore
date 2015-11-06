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

#ifndef __H__UG_refinement_mark_util_impl
#define __H__UG_refinement_mark_util_impl

#include <algorithm>
#include <limits>
#include "refinement_mark_util.h"
#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/refinement/refiner_interface.h"

namespace ug{

template <class TRef, class TIter, class TAAPos>
void MarkForAnisotropicRefinement (
			Grid& grid,
			TRef& ref,
			number minEdgeRatio,
			TIter elemsBegin,
			TIter elemsEnd,
			TAAPos aaPos)
{
	using namespace std;
	typedef typename TIter::value_type	elem_ptr_t;

	number minEdgeRatioSq = sq(minEdgeRatio);
	Grid::edge_traits::secure_container	assEdges;

	for(TIter iter = elemsBegin; iter != elemsEnd; ++iter){
		elem_ptr_t elem = *iter;
		ref.mark(elem, RM_CLOSURE);

		grid.associated_elements(assEdges, elem);

		if(assEdges.size() < 2)
			continue;

	//	find the length of the shortest and the longest edge
		number shortestLenSq = numeric_limits<number>::max();
		number longestLenSq = 0;

		for_each_in_vec(Edge* e, assEdges){
			number lenSq = EdgeLengthSq(e, aaPos);
			shortestLenSq = min(shortestLenSq, lenSq);
			longestLenSq = max(longestLenSq, lenSq);
		}end_for;

		if(longestLenSq < SMALL_SQ)
			continue;

		if(shortestLenSq / longestLenSq >= minEdgeRatioSq)
			continue;

	//	the element is anisotropic mark it and all long edges
		ref.mark(elem, RM_ANISOTROPIC);
	//	refine all edges that are at least half as long as the longest one
		number thresholdLenSq = shortestLenSq / minEdgeRatioSq;
		for_each_in_vec(Edge* e, assEdges){
			if(EdgeLengthSq(e, aaPos) > thresholdLenSq){
				ref.mark(e, RM_REFINE);
			}
		}end_for;
	}
}

}//	end of namespace

#endif	//__H__UG_refinement_mark_util_impl

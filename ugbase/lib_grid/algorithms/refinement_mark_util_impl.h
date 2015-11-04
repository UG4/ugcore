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

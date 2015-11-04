#ifndef __H__UG_refinement_mark_util
#define __H__UG_refinement_mark_util

namespace ug{

/**
 * \param grid	The whose elements shall be marked.
 *
 * \param ref	The specified refiner has to feature methods 'mark(elem_t*, byte)'
 *				where elem_t = Vertex, Edge, Face, and Volume. Note that you may also
 *				pass a selector instead of a refiner.
 *
 * \param minEdgeRatio	If the ratio between the shortest and the longest edge
 *						of an element is smaller than minEdgeRatio, the element
 *						is considered to be anisotropic. The element itself and
 *						all of its edges with a smaller ratio ar marked for refinement.
 *
 * \param elemsBegin	Iterator to the first element in the sequence of elements
 *						that shall be checked.
 *
 * \param elemsEnd		Iterator to the (pseudo-) element directly behind the
 *						last element in the sequence of elements that shall be checked.
 *
 * \param aaPos			A VertexAttachmentAccessor to an APosition compatible type.
 */
template <class TRef, class TIter, class TAAPos>
void MarkForAnisotropicRefinement (
			Grid& grid,
			TRef& ref,
			number minEdgeRatio,
			TIter elemsBegin,
			TIter elemsEnd,
			TAAPos aaPos);

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "refinement_mark_util_impl.h"

#endif	//__H__UG_refinement_mark_util

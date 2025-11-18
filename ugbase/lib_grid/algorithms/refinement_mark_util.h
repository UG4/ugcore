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


template <class TRef, class TEdgeIter, class TAAPos>
void MarkForRefinementByDirection (
			TRef& ref,
			TAAPos aaPos,
			TEdgeIter edgesBegin,
			TEdgeIter edgesEnd,
			const typename TAAPos::ValueType& dir,
			number minDeviationAngle,
			number maxDeviationAngle,
			bool selectFlipped);

}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "refinement_mark_util_impl.h"

#endif
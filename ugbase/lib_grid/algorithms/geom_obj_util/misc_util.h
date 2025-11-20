/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__MISC_UTIL__
#define __H__LIB_GRID__MISC_UTIL__

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid_objects/grid_objects.h"
#include "lib_grid/common_attachments.h"

namespace ug
{
/**
 * \brief contains miscellaneus methods that are related to GridObjects.
 * 
 * \defgroup lib_grid_algorithms_geom_obj_misc misc geometric object util
 * \ingroup lib_grid_algorithms
 * @{
 */


////////////////////////////////////////////////////////////////////////////////////////////
///	returns the shared side between the two elements or nullptr if no such side exists.
template <typename TElem>
UG_API
typename TElem::side*
GetSharedSide(Grid& grid, TElem* e1, TElem* e2);


////////////////////////////////////////////////////////////////////////
//	CalculateGridObjectCenter
///	calculates the center for arbitrary geometric object
template <typename TAAPosVRT>
UG_API
inline
typename TAAPosVRT::ValueType
CalculateGridObjectCenter(const GridObject* o, TAAPosVRT& aaPosVRT);


////////////////////////////////////////////////////////////////////////
///	returns the weighted center of the vertices of the given geometric object
/** TAAWeightVRT has to be an attachment to the vertices of the grid in which
 * the object is contained, with ValueType number (or compatible).
 */
template <typename TAAPosVRT, typename TAAWeightVRT>
UG_API
inline
typename TAAPosVRT::ValueType
CalculateGridObjectCenter(const GridObject* o, TAAPosVRT& aaPosVRT,
								TAAWeightVRT& aaWeight);

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
///	calculates the center for a set of elements
/**	TIterator::value_type has to be compatible with
 *  Vertex*, Edge*, Face* or Volume*.
 */
template <typename TIterator, typename TAAPosVRT>
UG_API 
typename TAAPosVRT::ValueType
CalculateCenter(TIterator begin, TIterator end, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	FindClosestByCoordinate
///	returns the element of the given type, whose center is closest to coord.
/**
 * Returns the element, whose center is closest to the given coordinate.
 *
 * TVertexPositionAttachmentAccessor has to be an AttachmentAccessor,
 * where AttachmentAccessor::ValueType is a vector-type compatible to
 * the lgmath vector descriptor.
 * The Accessor has to access an attachment of the vertices,
 * to which the faces between iterBegin and iterEnd refer.
 */
template<typename TElem, typename TVertexPositionAttachmentAccessor>
UG_API 
TElem* FindClosestByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
							   typename geometry_traits<TElem>::iterator iterBegin,
							   typename geometry_traits<TElem>::iterator iterEnd,
							   TVertexPositionAttachmentAccessor& aaPosVRT);

////////////////////////////////////////////////////////////////////////
///	Calculates the bounding box of a set of geometric objects
/**	TIterator has to be an iterator to a set containing elements of type
 * Edge*, Face* or Volume*. An overload for Vertex* exists.
 *
 * Make sure that TAAPos::ValueType == vector_t.
 */
template<typename vector_t, typename TIterator, typename TAAPos>
UG_API 
void CalculateBoundingBox(vector_t& vMinOut, vector_t& vMaxOut,
						  TIterator begin, TIterator end,
						  TAAPos& aaPos);

////////////////////////////////////////////////////////////////////////
//	NumSharedVertices
///	returns the number of vertices that are shared by two objects
/**	This algorithm uses Grid::mark.
 *
 *	Valid types are Edge*, Face*, Volume* and derivates of those.
 *	You may combine different types in one query.
 */
template <typename TElemPtr1, typename TElemPtr2>
UG_API 
size_t NumSharedVertices(Grid& grid, TElemPtr1 elem1, TElemPtr2 elem2);

////////////////////////////////////////////////////////////////////////
//	EraseConnectingElements
///	erases all elements that connect v1 and v2
UG_API 
void EraseConnectingElements(Grid& grid, Vertex* v1, Vertex* v2);

////////////////////////////////////////////////////////////////////////
//	EraseElements
///	erases all elements between iterBegin and iterEnd.
template <typename TElem>
UG_API 
void EraseElements(Grid& grid, typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd);

////////////////////////////////////////////////////////////////////////
//	ElementDiameter
///	returns the maximal squared distance between to element vertices
template <typename TElem, typename TAAPos>
number ElementDiameterSq(Grid& grid,
                         TAAPos& aaPos,
					     TElem* elem);

///	returns the maximal squared distance between to element vertices
template <typename TAAPos>
number ElementDiameterSq(Grid& grid,
                         TAAPos& aaPos,
					     GridObject* elem);

///	returns the maximal distance between to element vertices
template <typename TElem, typename TAAPos>
number ElementDiameter(Grid& grid,
                       TAAPos& aaPos,
					   TElem* elem);

///	returns the maximal diameter of all elements between iterBegin and iterEnd.
/** In parallel, the global max diameter is returned.*/
template <typename TAAPos, typename TIterator>
number MaxElementDiameter(Grid& grid, TAAPos& aaPos,
                          TIterator iterBegin, TIterator iterEnd);

///	returns the minimal diameter of all elements between iterBegin and iterEnd.
/** In parallel, the global min diameter is returned.*/
template <typename TAAPos, typename TIterator>
number MinElementDiameter(Grid& grid, TAAPos& aaPos,
                          TIterator iterBegin, TIterator iterEnd);


///	Returns the direction from the center of e1 to the center of e2
template <typename TElem1, typename TElem2, typename TAAPos>
typename TAAPos::ValueType
GetDirection (TElem1* e1, TElem2* e2, const TAAPos& aaPos);

///	Checks whether the center of e2 can be reached from the center of e1 in the given direction
/**
 * \param minAngle	minimal allowed deviation angle from 'dir' in degrees. Normally set to 0.
 * \param maxAngle	maximal allowed deviation angle from 'dir' in degrees */
template <typename TElem1, typename TElem2, typename TAAPos>
bool CheckDirection (TElem1* e1,
                     TElem2* e2,
                     const TAAPos& aaPos,
                     const typename TAAPos::ValueType& dir,
                     number minAngle,
                     number maxAngle);


/// @}				
}//	end of namespace

////////////////////////////////////////////////////////////////////////
//	include implementation of template methods.
#include "misc_util_impl.hpp"

#endif

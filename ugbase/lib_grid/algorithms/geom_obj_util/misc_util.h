//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d03

#ifndef __H__LIB_GRID__MISC_UTIL__
#define __H__LIB_GRID__MISC_UTIL__

#include "lib_grid/lg_base.h"

namespace ug
{
/**
 * \brief contains miscellaneus methods that are related to GeometricObjects.
 * 
 * \defgroup lib_grid_algorithms_geom_obj_misc misc geometric object util
 * \ingroup lib_grid_algorithms
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	CalculateCenter
///	calculates the center for a set of elements
/**	TIterator::value_type has to be compatible with
 *  VertexBase*, EdgeBase*, Face* or Volume*.
 */
template <class TIterator, class TAAPosVRT>
typename TAAPosVRT::ValueType
CalculateCenter(TIterator begin, TIterator end, TAAPosVRT& aaPos);

////////////////////////////////////////////////////////////////////////
//	FindByCoordinate
///	returns the element of the given type whose center is closest to coord.
/**
 * This method does not necessarily return the element that contains the given coordinate.
 * Instead it will simply search for the face whose center is closest to the specified
 * coordinate.
 * TVertexPositionAttachmentAccessor has to be an AttachmentAccessor,
 * where AttachmentAccessor::ValueType is a vector-type compatible to
 * the lgmath vector descriptor.
 * The Accessor has to access an attachment of the vertices,
 * to which the elements between iterBegin and iterEnd refer.
 *
 * This method can be called like this: FindByCoordinate<Face>(...).
 */
template<class TElem, class TVertexPositionAttachmentAccessor>
TElem* FindByCoordinate(const typename TVertexPositionAttachmentAccessor::ValueType& coord,
						typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd,
						TVertexPositionAttachmentAccessor& aaPosVRT);

////////////////////////////////////////////////////////////////////////
//	NumSharedVertices
///	returns the number of vertices that are shared by two objects
/**	This algorithm uses Grid::mark.
 *
 *	Valid types are EdgeBase*, Face*, Volume* and derivates of those.
 *	You may combine different types in one query.
 */
template <class TElemPtr1, class TElemPtr2>
size_t NumSharedVertices(Grid& grid, TElemPtr1 elem1, TElemPtr2 elem2);

////////////////////////////////////////////////////////////////////////
//	EraseConnectingElements
///	erases all elements that connect v1 and v2
void EraseConnectingElements(Grid& grid, VertexBase* v1, VertexBase* v2);

////////////////////////////////////////////////////////////////////////
//	EraseElements
///	erases all elements between iterBegin and iterEnd.
template <class TElem>
void EraseElements(Grid& grid, typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd);

/// @}				
}//	end of namespace

////////////////////////////////////////////////////////////////////////
//	include implementation of template methods.
#include "misc_util_impl.hpp"

#endif

//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d03

#ifndef __H__LIB_GRID__MISC_UTIL__
#define __H__LIB_GRID__MISC_UTIL__

#include "lib_grid/grid/grid.h"

namespace ug
{
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

////////////////////////////////////////////////////////////////////////
//	AssignIndices
///	assigns indices starting from 0.
/** Make sure that the accessor is valid and that it references an
 * attachment of the correct type. The elements between iterBegin and
 * iterEnd have to belong to the same grid for which the accessor was
 * created.*/
template <class TElem>
void AssignIndices(typename geometry_traits<TElem>::iterator iterBegin,
					typename geometry_traits<TElem>::iterator iterEnd,
					Grid::AttachmentAccessor<TElem, AInt>& aaInt);
				
}//	end of namespace

////////////////////////////////////////////////////////////////////////
//	include implementation of template methods.
#include "misc_util_impl.hpp"

#endif

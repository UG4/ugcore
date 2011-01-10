//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m07 d21

#ifndef __H__LIB_GRID__ATTACHMENT_UTIL__
#define __H__LIB_GRID__ATTACHMENT_UTIL__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/**
 * Several methods that ease attachment-handling are grouped here.
 *
 * \defgroup lib_grid_algorithms_attachment_util attachment util
 * \ingroup lib_grid_algorithms
 * @{
 */

////////////////////////////////////////////////////////////////////////
//	SetAttachmentValues
///	sets attachment-values for elements between elemsBegin and elemsEnd.
template <class TAttachmentAccessor, class TIter, class TVal>
void SetAttachmentValues(TAttachmentAccessor& aaVal,
						TIter elemsBegin, TIter elemsEnd,
						const TVal& val);

////////////////////////////////////////////////////////////////////////
//	ConvertMathVectorAttachmentValues
///	Fills the dest-attachment with values from the source-attachment.
/**
 * TSrcAttachment and TDestAttachment have to have ValueTypes that
 * are compatible with ug::MathVector.
 *
 * Copies values from srcAttachment to destAttachment.
 * if destAttachment is not already attached, it will be attached
 * automatically. The srcAttachment however has to be attached.
 *
 * Valid types for TElem are: VertexBase, EdgeBase, Face, Volume
 *
 * If the dimensions do not match, the algorithm behaves as follows:
 * dim(src) > dim(dest): Only dim(dest) values are copied per element.
 * dim(src) < dim(dest): Values in dimensions >= dim(src) are set to 0.
 */
template<class TElem, class TSrcAttachment, class TDestAttachment>
bool ConvertMathVectorAttachmentValues(Grid& grid,
							TSrcAttachment& srcAttachment,
							TDestAttachment& destAttachment);


////////////////////////////////////////////////////////////////////////
///	copies attachments from one grid to the other
/**
 * If aSrc is not attached to srcGrid, false is returned.
 * If aDest is not attached to destGrid, it is attached automatically.
 *
 * The method iterates through the elements specified by TElem
 * and copies the attachments.
 *
 * Call like this: CopyAttachments<VertexBase>(...);
 */
template <class TElem, class TAttachment>
bool CopyAttachments(Grid& srcGrid, TAttachment& aSrc,
					Grid& destGrid, TAttachment& aDest);

////////////////////////////////////////////////////////////////////////
///	copies attachments for the specified elements
/**
 * If aSrc is not attached to srcGrid, false is returned.
 * If aDest is not attached to destGrid, it is attached automatically.
 */
template <class TElemIter, class TAttachment>
bool CopyAttachments(Grid& grid, TElemIter elemsBegin, TElemIter elemsEnd,
					 TAttachment& aSrc, TAttachment& aDest);

////////////////////////////////////////////////////////////////////////
///	assigns indices to the elements between begin and end.
/**	Indices are stored in the given attachment. Make sure that the
 *	given attachment-accessor operates on an attachment-pipe at which
 *	those elements are registered.
 */
template <class TIterator, class TAAInt>
void AssignIndices(TIterator begin, TIterator end,
					TAAInt& aaInt, int baseIndex = 0);
							
/**@}*/ // end of doxygen defgroup command
}//	end of namespace

////////////////////////////////////////////////
// include implementations of template methods.
#include "attachment_util_impl.hpp"

#endif

//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d03

#ifndef __H__LIB_GRID__MISC_UTIL__IMPL__
#define __H__LIB_GRID__MISC_UTIL__IMPL__

#include "misc_util.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	EraseElements
template <class TElem>
void EraseElements(Grid& grid, typename geometry_traits<TElem>::iterator iterBegin,
						typename geometry_traits<TElem>::iterator iterEnd)
{
//	be careful to not invalidate the iterators.
	while(iterBegin != iterEnd)
	{
		TElem* e = *iterBegin;
		iterBegin++;
		grid.erase(e);
	}
}

}//	end of namespace

#endif

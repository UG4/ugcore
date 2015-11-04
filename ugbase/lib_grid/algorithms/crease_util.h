// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_crease_util
#define __H__UG_crease_util

#include "lib_grid/lg_base.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug{
template <class TVrtIter, class TAAPos>
void SelectKinkVertices(Grid& grid, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
						number thresholdAngle, bool selectDarts, TAAPos aaPos,
						Grid::edge_traits::callback cbConsiderEdge = ConsiderAll());
}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "crease_util_impl.h"

#endif	//__H__UG_crease_util

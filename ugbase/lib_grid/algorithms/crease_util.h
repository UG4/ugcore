// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_crease_util
#define __H__UG_crease_util

namespace ug{
template <class TVrtIter, class TAAPos>
void SelectKinkVertices(Grid& grid, TVrtIter vrtsBegin, TVrtIter vrtsEnd,
						number thresholdAngle, bool selectDarts, TAAPos aaPos,
						Grid::edge_traits::callback cbConsiderEdge =
							Grid::edge_traits::cb_consider_all);
}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "crease_util_impl.h"

#endif	//__H__UG_crease_util

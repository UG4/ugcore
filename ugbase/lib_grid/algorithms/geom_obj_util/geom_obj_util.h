//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m02 d02

#ifndef __H__LIB_GRID__GEOM_OBJ_UTIL__
#define __H__LIB_GRID__GEOM_OBJ_UTIL__

#include "vertex_util.h"
#include "edge_util.h"
#include "face_util.h"
#include "volume_util.h"
#include "misc_util.h"
#include "lib_grid/grid/grid_util.h"

namespace ug
{
///	calculates the center for a set of elements
/**	TIterator::value_type has to be compatible with
 *  VertexBase*, EdgeBase*, Face* or Volume*.
 */
template <class TIterator, class TAAPosVRT>
typename TAAPosVRT::ValueType
CalculateCenter(TIterator begin, TIterator end, TAAPosVRT& aaPos)
{
	int counter = 0;
	typename TAAPosVRT::ValueType center;
	VecSet(center, 0);
	for(TIterator iter = begin; iter != end; ++iter, ++counter)
		VecAdd(center, center, CalculateCenter(*iter, aaPos));
		
	if(counter > 0)
		VecScale(center, center, 1./(number)counter);
		
	return center;
}
}//	end of namespace

#endif

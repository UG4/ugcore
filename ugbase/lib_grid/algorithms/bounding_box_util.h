#ifndef __H__UG_bounding_box_util
#define __H__UG_bounding_box_util

#include "common/math/misc/shapes.h"

namespace ug{

///	calculates the smallest axis aligned box that contains the given vertex
template <class TAAPos>
AABox<typename TAAPos::ValueType>
CalculateBoundingBox (Vertex* e, TAAPos aaPos)
{
	return AABox<typename TAAPos::ValueType>(aaPos[e], aaPos[e]);
}

///	calculates the smallest axis aligned box that contains the given element
template <class TElem, class TAAPos>
AABox<typename TAAPos::ValueType>
CalculateBoundingBox (TElem* e, TAAPos aaPos)
{
	typedef AABox<typename TAAPos::ValueType> box_t;
	typename TElem::ConstVertexArray vrts = e->vertices();
	box_t box(aaPos[vrts[0]], aaPos[vrts[0]]);
	for(size_t i = 1; i < e->num_vertices(); ++i){
		box = box_t(box, aaPos[vrts[i]]);
	}
	return box;
}

}//	end of namespace

#endif	//__H__UG_bounding_box_util

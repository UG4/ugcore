// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 01.02.2012 (m,d,y)

#ifndef __H__UG__debug_util_impl__
#define __H__UG__debug_util_impl__

#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug
{

template <class TElem>
vector3 GetGeometricObjectCenter(Grid& g, TElem* elem)
{
	if(g.has_vertex_attachment(aPosition)){
		Grid::VertexAttachmentAccessor<APosition> aaPos(g, aPosition);
		return CalculateCenter(elem, aaPos);
	}
	else if(g.has_vertex_attachment(aPosition2)){
		Grid::VertexAttachmentAccessor<APosition2> aaPos(g, aPosition2);
		vector2 v = CalculateCenter(elem, aaPos);
		return vector3(v.x, v.y, 0);
	}
	if(g.has_vertex_attachment(aPosition1)){
		Grid::VertexAttachmentAccessor<APosition1> aaPos(g, aPosition1);
		vector1 v = CalculateCenter(elem, aaPos);
		return vector3(v.x, 0, 0);
	}

	UG_LOG("GetGeometricObjectCenter failed! No standard position attachment found.\n");
	return vector3(0, 0, 0);
}


template <class TElem>
int GetGeometricObjectIndex(Grid& g, TElem* elem)
{
	typedef typename Grid::traits<TElem>::base_object TBase;

	int counter = 0;
	for(typename Grid::traits<TBase>::iterator iter = g.begin<TBase>();
		iter != g.end<TBase>(); ++iter, ++counter)
	{
		if(*iter == elem)
			return counter;
	}
	return -1;
}

}//	end of namespace

#endif

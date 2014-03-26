// created by Sebastian Reiter
// s.b.reiter@gmail.com
// march 2014

#ifndef __H__UG_NORMAL_CALCULATION_IMPL__
#define __H__UG_NORMAL_CALCULATION_IMPL__

#include "normal_calculation.h"
#include "geom_obj_util/face_util.h"

namespace ug{

template <class TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Vertex* v, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType n;
	VecSet(n, 0);
	return n;
}

template <class TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Edge* e, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType n;
	VecSubtract(n, aaPos[e->vertex(sideIndex)],
				aaPos[e->vertex((sideIndex + 1) % 2)]);
	VecNormalize(n, n);
	return n;
}

template <class TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Face* f, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType c, n;
	c = CalculateCenter(f, aaPos);
	EdgeDescriptor ed;
	f->edge_desc(sideIndex, ed);

	DropAPerpendicular(n, c, aaPos[ed.vertex(0)], aaPos[ed.vertex(1)]);
	VecSubtract(n, n, c);
	VecNormalize(n, n);
	return n;
}

template <class TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Volume* v, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType n;
	FaceDescriptor fd;
	v->face_desc(sideIndex, fd);
	CalculateNormal(n, &fd, aaPos);
	return n;
}

template <class TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(GridObject* o, int sideIndex, TAAPos aaPos)
{
	int baseObjId = o->base_object_id();
	switch(baseObjId){
		case VERTEX:
			return CalculateOuterNormal(static_cast<Vertex*>(o), sideIndex, aaPos);
		case EDGE:
			return CalculateOuterNormal(static_cast<Edge*>(o), sideIndex, aaPos);
		case FACE:
			return CalculateOuterNormal(static_cast<Face*>(o), sideIndex, aaPos);
		case VOLUME:
			return CalculateOuterNormal(static_cast<Volume*>(o), sideIndex, aaPos);
	}
	UG_THROW("Unsupported base object id in CalculateOuterNormal");
}

}//	end of namespace

#endif

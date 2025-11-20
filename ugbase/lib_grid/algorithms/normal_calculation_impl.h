/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_NORMAL_CALCULATION_IMPL__
#define __H__UG_NORMAL_CALCULATION_IMPL__

#include "normal_calculation.h"
#include "geom_obj_util/face_util.h"
#include "common/math/ugmath.h"

namespace ug{

template <typename TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Vertex* v, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType n;
	VecSet(n, 0);
	return n;
}

template <typename TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Edge* e, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType n;
	VecSubtract(n, aaPos[e->vertex(sideIndex)],
				aaPos[e->vertex((sideIndex + 1) % 2)]);
	VecNormalize(n, n);
	return n;
}

template <typename TAAPos>
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

template <typename TAAPos>
inline typename TAAPos::ValueType
CalculateOuterNormal(Volume* v, int sideIndex, TAAPos aaPos)
{
	typename TAAPos::ValueType n;
	FaceDescriptor fd;
	v->face_desc(sideIndex, fd);
	CalculateNormal(n, &fd, aaPos);
	return n;
}

template <typename TAAPos>
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


inline vector2
CalculateNormal(EdgeVertices* edge,
				Grid::AttachmentAccessor<Vertex, Attachment<vector2> >& aaPos)
{
	vector2 d;
	VecSubtract(d, aaPos[edge->vertex(1)], aaPos[edge->vertex(0)]);
	VecNormalize(d, d);
	return vector2(d.y(), -d.x());
}


inline vector3
CalculateNormal(FaceVertices* face,
				Grid::AttachmentAccessor<Vertex, Attachment<vector3> >& aaPos)
{
	vector3 n;
	CalculateNormal(n, face, aaPos);
	return n;
}

}//	end of namespace

#endif

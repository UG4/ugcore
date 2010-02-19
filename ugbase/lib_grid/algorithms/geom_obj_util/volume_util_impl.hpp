//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m02 d19

#ifndef __H__LIB_GRID__VOLUME_UTIL_IMPL__
#define __H__LIB_GRID__VOLUME_UTIL_IMPL__

#include "lib_grid/lg_base.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron
inline bool
PointIsInsideTetrahedron(const vector3& v, Tetrahedron* tet,
						 Grid::VertexAttachmentAccessor<APosition>& aaPos)
{
	return PointIsInsideTetrahedron(v, aaPos[tet->vertex(0)], aaPos[tet->vertex(1)],
									aaPos[tet->vertex(2)], aaPos[tet->vertex(3)]);
}

}//	end of namespace

#endif

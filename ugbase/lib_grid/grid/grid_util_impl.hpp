//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

#ifndef __H__LIB_GRID__GRID_UTIL_IMPL__
#define __H__LIB_GRID__GRID_UTIL_IMPL__

#include <vector>
#include "grid_util.h"
#include "common/assert.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	EdgeContains
inline bool EdgeContains(EdgeVertices* e, VertexBase* vrt)
{
	return e->vertex(0) == vrt || e->vertex(1) == vrt;
}

inline bool EdgeContains(EdgeVertices* e, VertexBase* vrt1, VertexBase* vrt2)
{
	return ((e->vertex(0) == vrt1 && e->vertex(1) == vrt2)
			|| (e->vertex(1) == vrt1 && e->vertex(0) == vrt2));
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	new methods

////////////////////////////////////////////////////////////////////////
template <class TVrtContainer1, class TVrtContainer2>
bool CompareVertexContainer(const TVrtContainer1& con1,
					const TVrtContainer2& con2)
{
	int con1Size = (int)con1.size();

	if(con1Size != con2.size())
		return false;

	for(int i = 0; i < con1Size; ++i)
	{
		int j;
		for(j = 0; j < con1Size; ++j)
		{
			if(con1[i] == con2[j])
				break;
		}

	//	check whether we found a matching vertex
		if(j == con1Size)
			return false;
	}

	return true;
}


////////////////////////////////////////////////////////////////////////
//	GetVertex
inline VertexBase* GetVertex(VertexBase* vrt, size_t i)
{
	UG_ASSERT(i < 1, "A Vertex has only one vertex");
	return vrt;
}

inline VertexBase* GetVertex(EdgeBase* edge, size_t i)
{
	UG_ASSERT(i < edge->num_vertices(), "Wrong number of vertex");
	return edge->vertex(i);
}

inline VertexBase* GetVertex(Face* face, size_t i)
{
	UG_ASSERT(i < face->num_vertices(), "Wrong number of vertex");
	return face->vertex(i);
}

inline VertexBase* GetVertex(Volume* vol, size_t i)
{
	UG_ASSERT(i < vol->num_vertices(), "Wrong number of vertex");
	return vol->vertex(i);
}

}//	end of namespace libGrid

#endif

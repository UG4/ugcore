//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d16

#ifndef __H__LIB_GRID__GRID_UTIL_IMPL__
#define __H__LIB_GRID__GRID_UTIL_IMPL__

#include <vector>
#include "grid_util.h"
#include "common/assert.h"
#include "common/common.h"

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

inline size_t NumVertices(VertexBase* elem)
{
	return 1;
}

inline size_t NumVertices(EdgeBase* elem)
{
	return elem->num_vertices();
}

inline size_t NumVertices(Face* elem)
{
	return elem->num_vertices();
}

inline size_t NumVertices(Volume* elem)
{
	return elem->num_vertices();
}


////////////////////////////////////////////////////////////////////////
inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					  Grid& grid, VertexBase* v, bool clearContainer)
{
	CollectVertices(vVertexOut, grid, v, clearContainer);
}

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					   Grid& grid, EdgeBase* e, bool clearContainer)
{
	CollectVertices(vVertexOut, grid, e, clearContainer);
}

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					   Grid& grid, Face* f, bool clearContainer)
{
	CollectVertices(vVertexOut, grid, f, clearContainer);
}

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
					   Grid& grid, Volume* v, bool clearContainer)
{
	CollectVertices(vVertexOut, grid, v, clearContainer);
}

inline void CollectAssociated(std::vector<VertexBase*>& vVertexOut,
                              Grid& grid, GridObject* obj, bool clearContainer)
{
	uint type = obj->base_object_id();
	switch(type)
	{
		case VERTEX:CollectAssociated(vVertexOut, grid, reinterpret_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectAssociated(vVertexOut, grid, reinterpret_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectAssociated(vVertexOut, grid, reinterpret_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectAssociated(vVertexOut, grid, reinterpret_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}


////////////////////////////////////////////////////////////////////////
inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
					Grid& grid, VertexBase* vrt, bool clearContainer)
{
	CollectEdges(vEdgesOut, grid, vrt, clearContainer);
}

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
					Grid& grid, EdgeBase* e, bool clearContainer)
{
	CollectEdges(vEdgesOut, grid, e, clearContainer);
}

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
					Grid& grid, Face* f, bool clearContainer)
{
	CollectEdges(vEdgesOut, grid, f, clearContainer);
}

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
					Grid& grid, Volume* v, bool clearContainer)
{
	CollectEdges(vEdgesOut, grid, v, clearContainer);
}

inline void CollectAssociated(std::vector<EdgeBase*>& vEdgesOut,
                              Grid& grid, GridObject* obj, bool clearContainer)
{
	uint type = obj->base_object_id();
	switch(type)
	{
		case VERTEX:CollectAssociated(vEdgesOut, grid, reinterpret_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectAssociated(vEdgesOut, grid, reinterpret_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectAssociated(vEdgesOut, grid, reinterpret_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectAssociated(vEdgesOut, grid, reinterpret_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}


////////////////////////////////////////////////////////////////////////
inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, VertexBase* vrt, bool clearContainer)
{
	CollectFaces(vFacesOut, grid, vrt, clearContainer);
}

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, EdgeBase* e, bool clearContainer)
{
	CollectFaces(vFacesOut, grid, e, clearContainer);
}

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, Face* f, bool clearContainer)
{
	CollectFaces(vFacesOut, grid, f, clearContainer);
}

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
					Grid& grid, Volume* v, bool clearContainer)
{
	CollectFaces(vFacesOut, grid, v, clearContainer);
}

inline void CollectAssociated(std::vector<Face*>& vFacesOut,
                              Grid& grid, GridObject* obj, bool clearContainer)
{
	uint type = obj->base_object_id();
	switch(type)
	{
		case VERTEX:CollectAssociated(vFacesOut, grid, reinterpret_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectAssociated(vFacesOut, grid, reinterpret_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectAssociated(vFacesOut, grid, reinterpret_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectAssociated(vFacesOut, grid, reinterpret_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}


////////////////////////////////////////////////////////////////////////
inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, VertexBase* vrt, bool clearContainer)
{
	CollectVolumes(vVolumesOut, grid, vrt, clearContainer);
}

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, EdgeBase* e, bool clearContainer)
{
	CollectVolumes(vVolumesOut, grid, e, clearContainer);
}

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Face* f, bool clearContainer,
					bool ignoreAssociatedVolumes)
{
	CollectVolumes(vVolumesOut, grid, f, clearContainer, ignoreAssociatedVolumes);
}

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, Volume* vol, bool clearContainer)
{
	CollectVolumes(vVolumesOut, grid, vol, clearContainer);
}

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
					Grid& grid, FaceDescriptor& fd, bool clearContainer)
{
	CollectVolumes(vVolumesOut, grid, fd, clearContainer);
}

inline void CollectAssociated(std::vector<Volume*>& vVolumesOut,
                              Grid& grid, GridObject* obj, bool clearContainer)
{
	uint type = obj->base_object_id();
	switch(type)
	{
		case VERTEX:CollectAssociated(vVolumesOut, grid, reinterpret_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectAssociated(vVolumesOut, grid, reinterpret_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectAssociated(vVolumesOut, grid, reinterpret_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectAssociated(vVolumesOut, grid, reinterpret_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

inline void CollectVertices(std::vector<VertexBase*>& vVertexOut, Grid& grid,
                            				GridObject* obj, bool clearContainer)
{
	switch(obj->base_object_id())
	{
		case VERTEX:CollectVertices(vVertexOut, grid, static_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectVertices(vVertexOut, grid, static_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectVertices(vVertexOut, grid, static_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectVertices(vVertexOut, grid, static_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}

inline void CollectEdgesSorted(std::vector<EdgeBase*>& vEdgesOut, Grid& grid,
                            				GridObject* obj, bool clearContainer)
{
	switch(obj->base_object_id())
	{
		case VERTEX:CollectEdgesSorted(vEdgesOut, grid, static_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectEdgesSorted(vEdgesOut, grid, static_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectEdgesSorted(vEdgesOut, grid, static_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectEdgesSorted(vEdgesOut, grid, static_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}

inline void CollectFacesSorted(std::vector<Face*>& vFacesOut, Grid& grid,
                            				GridObject* obj, bool clearContainer)
{
	switch(obj->base_object_id())
	{
		case VERTEX:CollectFacesSorted(vFacesOut, grid, static_cast<VertexBase*>(obj), clearContainer); return;
		case EDGE:	CollectFacesSorted(vFacesOut, grid, static_cast<EdgeBase*>(obj), clearContainer); return;
		case FACE:	CollectFacesSorted(vFacesOut, grid, static_cast<Face*>(obj), clearContainer); return;
		case VOLUME:CollectFacesSorted(vFacesOut, grid, static_cast<Volume*>(obj), clearContainer); return;
	}
	throw(UGError("GeomObject type not known."));
}

}//	end of namespace libGrid

#endif

//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d19

#ifndef __H__LIBGRID__GRID_OBSERVERS__
#define __H__LIBGRID__GRID_OBSERVERS__

#include "common/types.h"

namespace ug
{
////////////////////////////////////////////////////////////////////////
//	predeclarations
class Grid;
class VertexBase;
class EdgeBase;
class Face;
class Volume;

////////////////////////////////////////////////////////////////////////
//	Observer types
enum ObserverType
{
	OT_NONE = 0,
	OT_GRID_OBSERVER = 1,
	OT_VERTEX_OBSERVER = 2,
	OT_EDGE_OBSERVER = 4,
	OT_FACE_OBSERVER = 8,
	OT_VOLUME_OBSERVER = 16,
	OT_FULL_OBSERVER = OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
						OT_FACE_OBSERVER | OT_VOLUME_OBSERVER
};

////////////////////////////////////////////////////////////////////////
//	GridObserver
class GridObserver
{
	public:
		virtual ~GridObserver()	{}

	//	grid callbacks
/*
		virtual void registered_at_grid(Grid* grid)			{}
		virtual void unregistered_from_grid(Grid* grid)		{}
*/
		virtual void grid_to_be_destroyed(Grid* grid)		{}
		virtual void elements_to_be_cleared(Grid* grid)		{}

	//	vertex callbacks
		virtual void vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent = NULL)	{}
		virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt)	{}
		virtual void vertex_to_be_replaced(Grid* grid, VertexBase* vrtOld, VertexBase* vrtNew)	{}

	//	edge callbacks
		virtual void edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent = NULL)		{}
		virtual void edge_to_be_erased(Grid* grid, EdgeBase* edge)		{}
		virtual void edge_to_be_replaced(Grid* grid, EdgeBase* edgeOld, EdgeBase* edgeNew)	{}

	//	face callbacks
		virtual void face_created(Grid* grid, Face* face, GeometricObject* pParent = NULL)				{}
		virtual void face_to_be_erased(Grid* grid, Face* face)			{}
		virtual void face_to_be_replaced(Grid* grid, Face* faceOld, Face* faceNew)	{}

	//	volume callbacks
		virtual void volume_created(Grid* grid, Volume* vol, GeometricObject* pParent = NULL)			{}
		virtual void volume_to_be_erased(Grid* grid, Volume* vol)		{}
		virtual void volume_to_be_replaced(Grid* grid, Volume* volOld, Volume* volNew)	{}
};

}//	end of namespace

#endif

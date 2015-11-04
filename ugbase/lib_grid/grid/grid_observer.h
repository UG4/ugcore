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
class Vertex;
class Edge;
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
/**
 * The grid observer defines an interface that can be specialized by
 * classes that want to be informed about changes in a grid.
 * If a class derives from GridObserver, it can be registered at a grid.
 * Registration is usually performed through a member function of the
 * observer class itself.
 * Most observers can only be registered at one grid at a time.
 *
 * Please note that methods of different observers are called in the
 * order in which they were registered at the grid. The only exception
 * are the vertex_to_be_erased, edge_to_be_erased, face_to_be_erased and
 * volume_to_be_erased. Those method are called in reverse order of
 * registration.
 */
class UG_API GridObserver
{
	public:
		virtual ~GridObserver()	{}

	//	grid callbacks
		virtual void grid_to_be_destroyed(Grid* grid)		{}
		virtual void elements_to_be_cleared(Grid* grid)		{}

	//	creation callbacks
	/**
	 *	\brief	Notified whenever a new element of the given type is created
	 *			in the given grid.
	 *
	 *	Creation callbacks are called in the order in which the GridObservers
	 * 	were registered at the given grid.
	 *
	 * 	If replacesParent is true, then pParent is of the same base type as the
	 * 	created object (e.g. in case of edge_created, the parent is an Edge*).
	 *  This case usually appears, when a contraining object is replaced by a
	 *  regular grid object if the same base type during refinement.
	 * 	The method is called with replacesParent == true by
	 * 	Grid::create_and_replace methods.
	 *
	 *	Please note: If replacesParent == true, then a call to
	 * 	OBJECT_to_be_erased(grid, pParent, obj) will follow (OBJECT
	 *  and obj are pseudonyms for the concrete type).*/
	/// \{
		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent = NULL,
									bool replacesParent = false)			{}

		virtual void edge_created(Grid* grid, Edge* e,
									GridObject* pParent = NULL,
									bool replacesParent = false)			{}

		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent = NULL,
									bool replacesParent = false)			{}

		virtual void volume_created(Grid* grid, Volume* vol,
									GridObject* pParent = NULL,
									bool replacesParent = false)			{}
	///	\}


	//	erase callbacks
	///	Notified whenever an element of the given type is erased from the given grid.
	/**	Erase callbacks are called in reverse order in which the GridObservers
	 * 	were registered at the given grid.
	 *
	 * 	if replacedBy != NULL the erased object is only replaced by another
	 *  grid object of the same base type. This usually happens when constraining
	 *  objects are replaced by regular objects in refinements. (E.g. a constraining
	 *  edge by become a regular Edge; note that both objects are of type
	 *  Edge*).
	 *
	 * \{ */
		virtual void vertex_to_be_erased(Grid* grid, Vertex* vrt,
										 Vertex* replacedBy = NULL)	{}

		virtual void edge_to_be_erased(Grid* grid, Edge* e,
										 Edge* replacedBy = NULL)	{}

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = NULL)	{}

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = NULL)	{}

	/**	\}	*/

	//	merge callbacks
	///	Notified when two elements of the same type are going to be merged.
	/**	Note that this method is invoked by Grid::objects_will_be_merged, which
	 * is called from outside the grid class. Implementors of algorithms in
	 * which objects are merged are thus responsible to call
	 * Grid::objects_will_be_merged.
	 *
	 * This callback is called in addition to ..._created and ..._to_be_erased
	 * callbacks and should thus only be used if small adjustments have to be
	 * made during a merge.
	 *
	 * Note that target may be identical to elem1 or elem2.
	 *
	 * \{ */
		virtual void vertices_to_be_merged(Grid* grid, Vertex* target,
										 Vertex* elem1, Vertex* elem2)	{}

		virtual void edges_to_be_merged(Grid* grid, Edge* target,
										 Edge* elem1, Edge* elem2)	{}

		virtual void faces_to_be_merged(Grid* grid, Face* target,
										 Face* elem1, Face* elem2)	{}

		virtual void volumes_to_be_merged(Grid* grid, Volume* target,
										 Volume* elem1, Volume* elem2)	{}

	/**	\}	*/
};

}//	end of namespace

#endif

/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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
enum ObserverType : uint8_t
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

using ObserverType_t = uint8_t;

constexpr ObserverType operator | (ObserverType lhs, ObserverType rhs) noexcept {
	return static_cast<ObserverType>(
		static_cast<ObserverType_t>(lhs) |
		static_cast<ObserverType_t>(rhs)
	);
}
constexpr ObserverType operator & (ObserverType lhs, ObserverType rhs) noexcept {
	return static_cast<ObserverType>(
		static_cast<ObserverType_t>(lhs) &
		static_cast<ObserverType_t>(rhs)
	);
}
constexpr ObserverType& operator |= (ObserverType &lhs, ObserverType rhs) noexcept {
		lhs = static_cast<ObserverType>(
			static_cast<ObserverType_t>(lhs) |
			static_cast<ObserverType_t>(rhs)
		);
	return lhs;
}


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
		virtual ~GridObserver()	= default;

	//	grid callbacks
		virtual void grid_to_be_destroyed(Grid* grid) {}
		virtual void elements_to_be_cleared(Grid* grid) {}

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
	 *  This case usually appears, when a containing object is replaced by a
	 *  regular grid object if the same base type during refinement.
	 * 	The method is called with replacesParent == true by
	 * 	Grid::create_and_replace methods.
	 *
	 *	Please note: If replacesParent == true, then a call to
	 * 	OBJECT_to_be_erased(grid, pParent, obj) will follow (OBJECT
	 *  and obj are pseudonyms for the concrete type).*/
	/// \{
		virtual void vertex_created(Grid* grid, Vertex* vrt,
									GridObject* pParent = nullptr,
									bool replacesParent = false) {}

		virtual void edge_created(Grid* grid, Edge* e,
									GridObject* pParent = nullptr,
									bool replacesParent = false) {}

		virtual void face_created(Grid* grid, Face* f,
									GridObject* pParent = nullptr,
									bool replacesParent = false) {}

		virtual void volume_created(Grid* grid, Volume* vol,
									GridObject* pParent = nullptr,
									bool replacesParent = false) {}
	///	\}


	//	erase callbacks
	///	Notified whenever an element of the given type is erased from the given grid.
	/**	Erase callbacks are called in reverse order in which the GridObservers
	 * 	were registered at the given grid.
	 *
	 * 	if replacedBy != nullptr the erased object is only replaced by another
	 *  grid object of the same base type. This usually happens when constraining
	 *  objects are replaced by regular objects in refinements. (E.g. a constraining
	 *  edge by become a regular Edge; note that both objects are of type
	 *  Edge*).
	 *
	 * \{ */
		virtual void vertex_to_be_erased(Grid* grid, Vertex* vrt,
										 Vertex* replacedBy = nullptr) {}

		virtual void edge_to_be_erased(Grid* grid, Edge* e,
										 Edge* replacedBy = nullptr) {}

		virtual void face_to_be_erased(Grid* grid, Face* f,
										 Face* replacedBy = nullptr) {}

		virtual void volume_to_be_erased(Grid* grid, Volume* vol,
										 Volume* replacedBy = nullptr) {}

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
										 Vertex* elem1, Vertex* elem2) {}

		virtual void edges_to_be_merged(Grid* grid, Edge* target,
										 Edge* elem1, Edge* elem2) {}

		virtual void faces_to_be_merged(Grid* grid, Face* target,
										 Face* elem1, Face* elem2) {}

		virtual void volumes_to_be_merged(Grid* grid, Volume* target,
										 Volume* elem1, Volume* elem2) {}

	/**	\}	*/
};

}//	end of namespace

#endif

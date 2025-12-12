/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_GIRD__DEBUG_UTIL__
#define __H__UG__LIB_GIRD__DEBUG_UTIL__

#include "lib_grid/lg_base.h"
#include "lib_grid/tools/surface_view.h"

namespace ug {

/**
 * Methods that log informations on grids, subset-handlers, ...
 * \defgroup lib_grid_algorithms_log_util log util
 * \ingroup lib_grid_algorithms
 * @{
 */

///	prints how many elements of each type exist in the goc.
void PrintElementNumbers(const GridObjectCollection& goc);

///	prints how many elements of each type exist in the grid.
void PrintGridElementNumbers(Grid& grid);

///	prints how many elements of each type exist in the multi grid.
void PrintGridElementNumbers(MultiGrid& mg);

///	prints how many elements of each type exist in the subset handler.
void PrintGridElementNumbers(GridSubsetHandler& sh);


///	prints information on all attachments of the specified grid
void PrintAttachmentInfo(Grid& grid);



///	Returns the center of the given element (SLOW - for debugging only!)
/**	Caution: This method is pretty slow and should only be used for debugging
 * purposes. It only works if one of the standard position attachments
 * aPosition, aPosition2 or aPosition1 is attached to the vertices of g.
 *
 * The method returns the center of the specified element in a vector3 structure,
 * regardless of the actual dimension of the position attachment. Unused
 * coordinates are set to 0.
 *
 * TElem has to be derived from GridObject
 *
 * \{
 */
template <typename TElem>
vector3 GetGridObjectCenter(Grid& g, TElem* elem);
inline vector3 GetGridObjectCenter(Grid& g, GridObject* elem);
/** \} */


/// Writes level and center of each object together with a custom value to a file
template <typename TElem, typename TAValue>
void WriteDebugValuesToFile(const char* filename, Grid& grid,
							TAValue& aVal, bool levelWise);


///	returns the index of the given element in the given grid.
/**	Runtime O(n). Returns -1 if the element could not be found.
 * \note: The all elements of the same base-type are considered during counting.*/
template <typename TElem>
int GetGridObjectIndex(Grid& g, TElem* elem);


///	Checks whether parent child connections in a multi-grid are correct
void CheckMultiGridConsistency(MultiGrid& mg);

///	checks whether all constraining and constrained objects are correctly connected.
bool CheckHangingNodeConsistency(Grid& g);

///	checks whether a multigrid is a valid haging vertex grid.
/**	This algorithm e.g. checks, whether adaptivity is represented by
 * constrained / constraining objects.
 *
 * Calls CheckHangingVertexConsistency(Grid&)
 */
bool CheckHangingNodeConsistency(MultiGrid& mg);



///	Checks whether distributed objects have the same type on all processes.
/**	Especially when constrained / constraining objects are present in a grid,
 * it may be useful to check the types if errors occur during refinement / distribution.
 *
 * \return	true if all object types match, false if not.
 */
bool CheckDistributedObjectConstraintTypes(MultiGrid& mg);


///	Check whether local parent types match the type of the actual parent element.
bool CheckDistributedParentTypes(MultiGrid& mg);


///	Checks whether associated elements and associated constrained/constraining objects are fine
/** \{ */
bool CheckElementConsistency(MultiGrid& mg, Vertex* v);
bool CheckElementConsistency(MultiGrid& mg, Edge* e);
bool CheckElementConsistency(MultiGrid& mg, Face* f);
/** \} */


///	Returns a string containing information on the given element
/**	The string contains the position and level of the element, whether it is
 * normal, constrained or constraining, and its assigned and actual parent types.
 * In parallel the interface states are also returned.
 *
 * This method is intended for use in UG_ASSERT or UG_THROW to print additional
 * information on an element for which a problem occurred.
 * \return	string containing gathered information on the given element
 * \{ */
std::string ElementDebugInfo(const Grid& grid, GridObject* e);
/** \} */

///	Performs some tests on a surface-view (checks iterators vs surface-states)
//void CheckSurfaceViewConsistency(SurfaceView& sv);

/**@}*/ // end of doxygen defgroup command
}//	end of namespace

////////////////////////////////////////
//	include implementation
#include "debug_util_impl.hpp"

#endif

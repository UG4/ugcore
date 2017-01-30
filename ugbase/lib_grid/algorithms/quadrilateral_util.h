/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG_quadrilateral_util
#define __H__UG_quadrilateral_util

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid_objects/grid_objects_2d.h"
#include "common/error.h"

namespace ug{

///	Creates a new quadrilateral from the specified triangles but doesn't register it at the grid.
/** Ownership of the returned quadrilateral is transferred to the caller.
 * The two triangles have to share exactly 2 vertices.*/
UG_API
Quadrilateral* CreateQuadrilateral_NoRegistration(Grid& g, Face* tri1, Face* tri2);

///	Creates a new quadrilateral from the specified triangles in the specified grid
/** The two triangles have to share exactly 2 vertices.*/
UG_API
Quadrilateral* CreateQuadrilateral(Grid& g, Face* tri1, Face* tri2);

///	Replaces the specified triangles by one quadrilateral.
/** The two triangles have to share exactly 2 vertices.*/
UG_API
Quadrilateral* ReplaceByQuadrilateral(Grid& g, Face* tri1, Face* tri2);

///	Attempts to replace the given set of triangles by a set of quadrilaterals
template <class face_iter_t, class TAAPos>
void ReplaceByQuadrilaterals_FaceBased(
		Grid& g,
		face_iter_t facesBegin,
		face_iter_t facesEnd,
		TAAPos aaPos);

///	Attempts to replace triangles associated with the given set of edges by a set of quadrilaterals
/**	Quadrilaterals are generated in the order from good quads to bad quads.*/
template <class edge_iter_t, class TAAPos>
void ReplaceByQuadrilaterals_EdgeBased(
		Grid& g,
		edge_iter_t edgesBegin,
		edge_iter_t edgesEnd,
		TAAPos aaPos);


///	Attempts to replace the given set of faces by a set of quadrilaterals
/** The method finds edges which connect two triangles in the specified sequence
 * and will pass those to 'ReplaceByQuadrilaterals_EdgeBasedNoSort'
 * Quadrilaterals are generated in the order from good quads to bad quads.*/
template <class face_iter_t>
void ReplaceByQuadrilaterals_FaceBasedNoSort(
		Grid& g,
		face_iter_t facesBegin,
		face_iter_t facesEnd);

///	Attempts to replace triangles associated with the given set of edges by a set of quadrilaterals
/** The method tries for each edge in the specified order to replace associated
 * triangles by one quadrilateral. The replace is only performed if exactly
 * two different triangles are connected with the edge. Otherwise the edge
 * is simply ignored.
 *
 * \warning	Each edge may only be contained once in the given sequence.*/
template <class edge_iter_t>
void ReplaceByQuadrilaterals_EdgeBasedNoSort(
		Grid& g,
		edge_iter_t edgesBegin,
		edge_iter_t edgesEnd);

}//	end of namespace


////////////////////////////////////////
//	include implementation
#include "quadrialteral_util_impl.h"

#endif	//__H__UG_quadrilateral_util

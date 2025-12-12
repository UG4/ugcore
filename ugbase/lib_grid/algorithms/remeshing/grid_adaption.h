/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_GRID__GRID_ADAPTION__
#define __H__LIB_GRID__GRID_ADAPTION__

// #include <vector>
#include "lib_grid/lg_base.h"

namespace ug {

///	\addtogroup lib_grid_algorithms_remeshing
///	@{

////////////////////////////////////////////////////////////////////////
///	Performs local remeshing so that the grid is adapted to the given cylinder.
/**
 * This algorithm uses Grid::mark.
 *
 * The resulting grid features an edge-loop that approximates the
 * intersection of the cylinder with the grid-surface.
 * New faces are inserted accordingly.
 *
 * Please note that volumes-geometries are not supported.
 * 
 * The given selector (sel) will contain the faces that lie inside
 * the cylinder when the algorithm is done.
 *
 * aInt has to be attached to the edges of the grid.
 * It will be used to store temporary values. Initial values are ignored.
 *
 * The algorithm requires the option FACEOPT_AUTOGENERATE_EDGES in grid.
 * If it isn't enabled it will be automatically enabled.
 *
 * An overloaded and somewhat slower version of this method exists,
 * which automatically creates the aInt attachment. It is thus
 * a little more comfortable for the caller.
 *
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 */
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
						   Vertex* vrtCenter, const vector3& normal,
						   number radius, number rimSnapThreshold, AInt& aInt,
						   APosition& aPos = aPosition);

////////////////////////////////////////////////////////////////////////
///	Performs local remeshing so that the grid is adapted to the given cylinder.
/**
 * This is an overloaded version of AdaptSurfaceGridToCylinder.
 * You don't have to pass an edge-integer-attachment to this version,
 * which makes it a little more comfortable, but at the same time
 * slower. You should worry about this slow-down if you intend to
 * call this method multiple times. You should then consider
 * to call the faster version of AdaptSurfaceGridToCylinder,
 * which takes aInt as a parameter.
 *
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 */
bool AdaptSurfaceGridToCylinder(Selector& selOut, Grid& grid,
						   Vertex* vrtCenter, const vector3& normal,
						   number radius, number rimSnapThreshold,
						   APosition& aPos = aPosition);

/// @}	// end of add_to_group command

}//	end of namespace

#endif

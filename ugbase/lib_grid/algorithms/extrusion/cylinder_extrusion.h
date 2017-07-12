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

#ifndef __H__LIB_GRID__CYLINDER_EXTRUSION__
#define __H__LIB_GRID__CYLINDER_EXTRUSION__

#include "lib_grid/lg_base.h"

namespace ug
{
/// \addtogroup lib_grid_algorithms_extrusion
///	@{

////////////////////////////////////////////////////////////////////////
///	adapts the grid around the given vertex to a cylinder and extrudes it.
/**
 * If you're not interested in the subsets of the extruded geometry, or if
 * you want to handle them yourself, please use a overloaded version of
 * this method.
 * \todo: 	Add a parameter that allows to specify whether the
 *	 		source-faces of the extrusion shall be erased.
 *
 * \param grid: The Grid
 * \param sh: 	The SubsetHandler
 * \param vrt: 	The Vertex
 * \param direction: the direction
 *
 * \param height: 	The actual extrude-amount is determined by scaling direction
 *					with height. Note that if height is negative, the orientation
 *					of the extruded geometrie will be inverted.
 *
 * \param radius: radius
 *
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 *
 * \param aaPos: position attachment accessor
 *
 * \param bottomSubInd: default value is -1. Defines the subset into which
 *						the bottom-faces of the cylinder will go
 *
 * \param cylSubInd: 	default value is -1. Defines the subset into which
 *						the bottom-faces of the cylinder will go
 *
 * \param pSel:	You may pass a pointer to a selector that may be used by
 *				this method internally. This makes sense if you call
 *				this method repeatedly, since a repeated allocation and
 *				deallocation can be avoided.
 */
bool ExtrudeCylinder(Grid& grid, SubsetHandler& sh, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd = -1, int cylSubInd = -1,
					Selector* pSel = NULL);

///	adapts the grid around the given vertex to a cylinder and extrudes it.
/**
 * If you want the method to automatically assign subset-indices to the
 * extruded geometry, then please use an overloaded version of this method.
 *
 * \todo: 	Add a parameter that allows to specify whether the
 *	 		source-faces of the extrusion shall be erased.
 *
 * \param grid: The Grid
 * \param vrt: 	The Vertex
 * \param direction: the direction
 * \param radius: radius
 * \param rimSnapThreshold:	If a vertex lies closer to the rim than rimSnapThreshold,
 * 							then it will be projected to the rim.
 * \param aaPos: position attachment accessor
 *
 * \param height: 	The actual extrude-amount is determined by scaling direction
 *					with height. Note that if height is negative, the orientation
 *					of the extruded geometrie will be inverted.
 *
 * \param pSel:	You may pass a pointer to a selector that may be used by
 *				this method internally. This makes sense if you call
 *				this method repeatedly, since a repeated allocation and
 *				deallocation can be avoided.
 *
 * \param minDot:	Faces whose normal has a dot-product lower that minDot
 *					whith the given direction are not regarded as cylinder-
 *					bottom faces and are thus not extruded.
 */					
bool ExtrudeCylinder(Grid& grid, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Selector* pSel = NULL);
					
/// @}

}//	end of namespace

#endif

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

#ifndef __H__LIB_GRID__REGULAR_REFINER__
#define __H__LIB_GRID__REGULAR_REFINER__

#include "lib_grid/lg_base.h"
#include "refinement_callbacks.h"
#include "common/ug_config.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_refinement
///	@{

////////////////////////////////////////////////////////////////////////
//	Refine
///	refines selected faces and edges regularily and builds a closure on adjacent unselected faces.
/**
 * Selected faces will be refined with regular refinement.
 * Faces which are adjacent to selected faces or selected edges and
 * which are not selected, will be refined too to build a closure.
 * In this case however, the refinement is not regular - faces are
 * refined in a way to avoid hanging nodes.
 *
 * Pass a grid and a selector (which is working on the grid). The
 * provided aInt is required by the algorithm to store temporary
 * values. You have to pass it to the algorithm to allow
 * maximal speed for repeated small refinements.
 * 
 * All involved geometric objects will be selected after the method
 * terminated (This includes vertices adjacent to selected edges and faces).
 *
 * aInt has to be attached to the edges of the grid.
 *
 * If you are interested in rare, big refinements, you may also use
 * the overloaded version of Refine, which only takes a Grid and a 
 * Selector.
 *
 * If 'useSnapPoints' is set to 'true' (default is 'false'), selected vertices
 * are used as snap-points. They play a role during refinement of quadrilaterals
 * and of volume-elements with quadrilateral sides.
 *
 * \sa ug::RegularRefiner, ug::HangingNodeRefiner
 */
UG_API
bool Refine(Grid& grid, Selector& sel, AInt& aInt,
			IRefinementCallback* refCallback = NULL,
			bool useSnapPoints = false);

///	refines selected faces and edges regularily and builds a closure on adjacent unselected faces.
/**
 * This method overloads Refine(Grid& grid, Selector& sel, AInt& aInt).
 * It is slower that the full version of refine, since it tempiroraily
 * attaches the required aInt attachment to the edges of aInt before
 * it calls the original version.
 *
 * This method should only be used if only very few refinement steps
 * are performed and if speed is not cruical.
 *
 * \sa ug::RegularRefiner, ug::HangingNodeRefiner
 */
UG_API
bool Refine(Grid& grid, Selector& sel,
			IRefinementCallback* refCallback = NULL,
			bool useSnapPoints = false);

/// @}
}// end of namespace

#endif

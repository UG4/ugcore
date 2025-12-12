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

#ifndef __H__UG__hanging_node_refiner_t__
#define __H__UG__hanging_node_refiner_t__

#include "hanging_node_refiner_grid.h"
#include "hanging_node_refiner_multi_grid.h"

namespace ug {


///	Gives access to a hanging node refiner, depending on the grid-type
template <typename TGrid>
class THangingNodeRefiner;

template <>
class THangingNodeRefiner<Grid> : public HangingNodeRefiner_Grid
{
	public:
	explicit THangingNodeRefiner(SPRefinementProjector projector = nullptr) :
			HangingNodeRefiner_Grid(projector)	{}

	explicit THangingNodeRefiner(Grid& grid, SPRefinementProjector projector = nullptr) :
			HangingNodeRefiner_Grid(grid, projector)	{}
};

template <>
class THangingNodeRefiner<MultiGrid> : public HangingNodeRefiner_MultiGrid
{
	public:
	explicit THangingNodeRefiner(SPRefinementProjector projector = nullptr) :
			HangingNodeRefiner_MultiGrid(projector)	{}

	explicit THangingNodeRefiner(MultiGrid& mg, SPRefinementProjector projector = nullptr) :
			HangingNodeRefiner_MultiGrid(mg, projector)	{}
};

}//	end of namespace

#endif

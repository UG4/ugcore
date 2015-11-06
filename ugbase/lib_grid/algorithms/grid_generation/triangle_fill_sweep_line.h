/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include <list>
#include <vector>
#include <cassert>
#include <map>
#include <stack>
#include "lib_grid/lg_base.h"
#include "common/ug_config.h"

#ifndef __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE__
#define __H__UG4__LIB_GRID__TRIANGLE_FILL_SWEEP_LINE__

namespace ug
{

///	Performs triangulation of a polygon and resoves additional inner edges.
/**	This algorithm triangulates a set of edges.
 *	It is important, that a closed outer boundary exists.
 *	All edges are resolved in the final geometry.*/
UG_API
bool TriangleFill_SweepLine(std::vector<int>& facesOut,
							const std::vector<vector2>& srcVrts,
							/*const */std::vector<int>& srcEdges);

///	Performs triangulation of a 3d polygon and resolves inner edges.
/**	The polygon should lie in a 2d hyperplane.*/
UG_API
bool TriangleFill_SweepLine(std::vector<int>& facesOut,
							const std::vector<vector3>& srcVrts,
							/*const */std::vector<int>& srcEdges);

///	Performs triangulation of a polygon and resoves additional inner edges.
/**
 *	This algortighm uses Grid::mark.
 *
 *	The polygon should lie in a 2d hyperplane.
 */
template <class TIterator>
bool TriangleFill_SweepLine(Grid& grid, TIterator edgesBegin,
							TIterator edgesEnd, APosition& aPosVRT,
							AInt& aIntVRT,
							SubsetHandler* pSH = NULL,
							int newSubsetIndex = -1);

}//	namespace ug

////////////////////////////////
//	include implementation
#include "triangle_fill_sweep_line_impl.hpp"

#endif

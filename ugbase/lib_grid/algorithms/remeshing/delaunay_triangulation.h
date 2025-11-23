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

#ifndef __H__UG__delaunay_triangulation__
#define __H__UG__delaunay_triangulation__

// #include <queue>
// #include <vector>
// #include <sstream>
#include "delaunay_info.h"
// #include "common/ug_config.h"
#include "lib_grid/lg_base.h"
// #include "lib_grid/algorithms/geom_obj_util/edge_util.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
template <typename TAAPos>
bool MakeDelaunay(DelaunayInfo<TAAPos>& info);

////////////////////////////////////////////////////////////////////////////////
///	Transforms the given triangle-set into a delaunay set
/** Creates a delaunay triangulation. If a minAngle greater than 0 is specified,
 * then additional vertices are introduced, if required to generate fulfill
 * the min-angle-condition.
 *
 * Make sure that init_marks was performed on the specified deaunay-info object.
 * If the triangulation is not already delaunay, all interior edges should also
 * be pushed to the candidates pool (can be done automatically in DelaunayInfo::init_marks)
 */
template <typename TAAPos>
bool QualityGridGeneration(Grid& grid, DelaunayInfo<TAAPos>& info,
						   number minAngle = 0,
				  	  	   int maxSteps = -1/*remove this*/);

template <typename TriIter, typename TAAPos>
bool QualityGridGeneration(Grid& grid, TriIter trisBegin, TriIter trisEnd,
						   TAAPos& aaPos, number minAngle = 0,
				  	  	   Grid::edge_traits::callback cbConstrainedEdge = ConsiderNone(),
				  	  	   int maxSteps = -1/*remove this*/)
{
//	set up a delaunay-info structure
	DelaunayInfo<TAAPos> info(grid, aaPos, cbConstrainedEdge);
	info.init_marks(trisBegin, trisEnd, true);
	return QualityGridGeneration(grid, info, minAngle, maxSteps);	
}

}//	end of namespace

#endif

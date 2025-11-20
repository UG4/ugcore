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

#ifndef __H__LIB_GRID__EXTRUDER_UTIL_IMPL__
#define __H__LIB_GRID__EXTRUDER_UTIL_IMPL__

#include <vector>
#include "extruder_util.h"
#include "extrusion/extrude.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
template <typename TIterator>
bool RepeatedVertexExtrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const ug::vector3& stepDir)
{
	std::vector<Vertex*> vElems;
	vElems.assign(iterBegin, iterEnd);

	for(int i = 0; i < numSteps; ++i)
		Extrude(grid, &vElems, nullptr, nullptr, stepDir);

	return true;
}

////////////////////////////////////////////////////////////////////////
template <typename TIterator>
bool RepeatedEdgeExtrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const ug::vector3& stepDir)
{
	std::vector<Edge*> vElems;
	vElems.assign(iterBegin, iterEnd);

	for(int i = 0; i < numSteps; ++i)
		Extrude(grid, nullptr, &vElems, nullptr, stepDir);

	return true;
}

////////////////////////////////////////////////////////////////////////
template <typename TIterator>
bool RepeatedFaceExtrusion(Grid& grid,
							TIterator iterBegin, TIterator iterEnd,
							int numSteps, const ug::vector3& stepDir)
{
	std::vector<Face*> vElems;
	vElems.assign(iterBegin, iterEnd);

	for(int i = 0; i < numSteps; ++i)
		Extrude(grid, nullptr, nullptr, &vElems, stepDir);

	return true;
}

}//	end of namespace

#endif

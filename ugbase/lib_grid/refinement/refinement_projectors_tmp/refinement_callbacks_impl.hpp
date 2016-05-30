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


#ifndef __H__LIB_GRID__REFINEMENT_CALLBACKS_IMPL__
#define __H__LIB_GRID__REFINEMENT_CALLBACKS_IMPL__

#include <cassert>
#include <algorithm>
#include "refinement_callbacks.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/subdivision/subdivision_rules_piecewise_loop.h"

namespace ug
{

template <class TAttachmentAccessor>
int IRefinementCallback::
current_pos_helper(number* coordsOut, Vertex* vrt, int maxCoords,
				   TAttachmentAccessor& aaPos)
{
	using namespace std;
	const int numCoords = min(maxCoords, (int)TAttachmentAccessor::ValueType::Size);
	typename TAttachmentAccessor::ValueType& p = aaPos[vrt];
	for(int i = 0; i < numCoords; ++i)
		coordsOut[i] = p[i];
	return numCoords;
}

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class TAPosition>
RefinementCallbackLinear<TAPosition>::
RefinementCallbackLinear() :
	m_pGrid(NULL)
{
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
RefinementCallbackLinear<TAPosition>::
RefinementCallbackLinear(Grid& grid, TAPosition& aPos) :
	m_pGrid(&grid)
{
//	we have to make sure that aPos is attached at the grid.
//	This is important to avoid crashes later on.
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	m_aaPos.access(grid, aPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
RefinementCallbackLinear<TAPosition>::
~RefinementCallbackLinear()
{
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = m_aaPos[parent];
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
int RefinementCallbackLinear<TAPosition>::
current_pos(number* coordsOut, Vertex* vrt, int maxCoords)
{
	return IRefinementCallback::current_pos_helper(coordsOut, vrt, maxCoords, m_aaPos);
}




}// end of namespace

#endif


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

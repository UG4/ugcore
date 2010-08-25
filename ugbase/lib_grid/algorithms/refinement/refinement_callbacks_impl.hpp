// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m08 d25


#ifndef __H__LIB_GRID__REFINEMENT_CALLBACKS_IMPL__
#define __H__LIB_GRID__REFINEMENT_CALLBACKS_IMPL__

#include <cassert>
#include "refinement_callbacks.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug
{

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
new_vertex(VertexBase* vrt, VertexBase* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = m_aaPos[parent];
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(VertexBase* vrt, EdgeBase* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(VertexBase* vrt, Face* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void RefinementCallbackLinear<TAPosition>::
new_vertex(VertexBase* vrt, Volume* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

}// end of namespace

#endif

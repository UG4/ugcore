// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m07 d08

#include <cassert>
#include "refinement_callbacks.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
RefinementCallbackLinear::
RefinementCallbackLinear() :
	m_pGrid(NULL)
{
}

////////////////////////////////////////////////////////////////////////
RefinementCallbackLinear::
RefinementCallbackLinear(Grid& grid, APosition& aPos) :
	m_pGrid(&grid)
{
	m_aaPos.access(grid, aPos);
}

////////////////////////////////////////////////////////////////////////
RefinementCallbackLinear::
~RefinementCallbackLinear()
{
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackLinear::
new_vertex(VertexBase* vrt, VertexBase* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = m_aaPos[parent];
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackLinear::
new_vertex(VertexBase* vrt, EdgeBase* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackLinear::
new_vertex(VertexBase* vrt, Face* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackLinear::
new_vertex(VertexBase* vrt, Volume* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = CalculateCenter(parent, m_aaPos);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
RefinementCallbackEdgePlaneCut::
RefinementCallbackEdgePlaneCut()
{
}

////////////////////////////////////////////////////////////////////////
RefinementCallbackEdgePlaneCut::
RefinementCallbackEdgePlaneCut(Grid& grid, const vector3& p,
										const vector3& n,
										APosition& aPos) :
	RefinementCallbackLinear(grid, aPos),
	m_p(p),
	m_n(n)
{
}

////////////////////////////////////////////////////////////////////////
RefinementCallbackEdgePlaneCut::
~RefinementCallbackEdgePlaneCut()
{
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackEdgePlaneCut::
new_vertex(VertexBase* vrt, EdgeBase* parent)
{
	number t;
	vector3 v;
	vector3 dir;
	VecSubtract(dir, m_aaPos[parent->vertex(1)], m_aaPos[parent->vertex(0)]);
	
	if(RayPlaneIntersection(v, t, m_aaPos[parent->vertex(0)], dir, m_p, m_n))
	{
		m_aaPos[vrt] = v;
	}
	else{
		RefinementCallbackLinear::new_vertex(vrt, parent);
	}
}

}// end of namespace

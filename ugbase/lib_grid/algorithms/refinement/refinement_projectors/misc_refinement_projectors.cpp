#include <cassert>
#include <algorithm>
#include "misc_refinement_projectors.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
using namespace std;

namespace ug
{

////////////////////////////////////////////////////////////////////////
RefinementCallback_IntersectCylinder::
RefinementCallback_IntersectCylinder()
{}

///	make sure that aPos is attached to the vertices of the grid.
RefinementCallback_IntersectCylinder::
RefinementCallback_IntersectCylinder(Grid& grid, const vector3& center,
									 const vector3& axis, number radius,
									 APosition& aPos) :
	 RefinementCallbackLinear<APosition>(grid, aPos),
	 m_center(center),
	 m_axis(axis),
	 m_radius(radius)
{}

RefinementCallback_IntersectCylinder::
~RefinementCallback_IntersectCylinder()
{}

void RefinementCallback_IntersectCylinder::
new_vertex(Vertex* vrt, Edge* parent)
{
	if(parent){
		number t0, t1;
		vector3 from = m_aaPos[parent->vertex(0)];
		vector3 dir;
		VecSubtract(dir, m_aaPos[parent->vertex(1)], from);

		if(RayCylinderIntersection(t0, t1, from, dir, m_center, m_axis, m_radius))
		{
		//	if there are two intersections with parameters between 0 and 1,
		//	we'll return their median.
			bool t0IsFine = (t0 >= 0) && (t0 <= 1);
			bool t1IsFine = (t1 >= 0) && (t1 <= 1);
			if(t0IsFine){
				if(t1IsFine){
					vector3 v0, v1;
					VecScaleAdd(v0, 1., from, t0, dir);
					VecScaleAdd(v1, 1., from, t1, dir);
					VecScaleAdd(m_aaPos[vrt], 0.5, v0, 0.5, v1);
				}
				else
					VecScaleAdd(m_aaPos[vrt], 1., from, t0, dir);
				return;
			}
			else if(t1IsFine){
				VecScaleAdd(m_aaPos[vrt], 1., from, t1, dir);
				return;
			}
		}
	}
	RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);
}


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
	RefinementCallbackLinear<APosition>(grid, aPos),
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
new_vertex(Vertex* vrt, Edge* parent)
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
		RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);
	}
}





}// end of namespace

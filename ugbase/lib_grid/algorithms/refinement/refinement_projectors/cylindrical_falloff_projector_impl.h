// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__cylindrical_falloff_projector_impl__
#define __H__UG__cylindrical_falloff_projector_impl__

#include "cylindrical_falloff_projector.h"

namespace ug{

template <class TAPosition>
CylindricalFalloffProjector<TAPosition>::
CylindricalFalloffProjector() :
	m_pGrid(NULL)
{
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
CylindricalFalloffProjector<TAPosition>::
CylindricalFalloffProjector(Grid& grid, TAPosition& aPos,
						 const typename TAPosition::ValueType& center,
						 const typename TAPosition::ValueType& axis,
						 number innerRadius, number outerRadius) :
	m_pGrid(&grid),
	m_center(center),
	m_axis(axis),
	m_innerRadius(innerRadius),
	m_outerRadius(outerRadius)
{
//	we have to make sure that aPos is attached at the grid.
//	This is important to avoid crashes later on.
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	m_aaPos.access(grid, aPos);

	if(m_outerRadius < m_innerRadius + SMALL)
		m_outerRadius = m_innerRadius + SMALL;
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
CylindricalFalloffProjector<TAPosition>::
~CylindricalFalloffProjector()
{
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CylindricalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = m_aaPos[parent];
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CylindricalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	perform_projection(vrt, parent);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CylindricalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
	perform_projection(vrt, parent);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
void CylindricalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
	perform_projection(vrt, parent);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
int CylindricalFalloffProjector<TAPosition>::
current_pos(number* coordsOut, Vertex* vrt, int maxCoords)
{
	return IRefinementCallback::current_pos_helper(coordsOut, vrt, maxCoords, m_aaPos);
}

////////////////////////////////////////////////////////////////////////
template <class TAPosition>
template <class TElem>
void CylindricalFalloffProjector<TAPosition>::
perform_projection(Vertex* vrt, TElem* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");

//	calculate the new position by linear interpolation and project that point
//	onto the cylinder.
	typename TElem::ConstVertexArray vrts = parent->vertices();
	size_t numVrts = parent->num_vertices();
	number avDist = 0;
	pos_type parentCenter;
	VecSet(parentCenter, 0);

	for(size_t i = 0; i < numVrts; ++i){
		const pos_type& p = m_aaPos[vrts[i]];
		avDist += DistancePointToRay(p, m_center, m_axis);
		VecAdd(parentCenter, parentCenter, p);
	}

	avDist /= (number)numVrts;
	VecScale(parentCenter, parentCenter, 1. / (number)numVrts);

	number ia = (avDist - m_innerRadius) / (m_outerRadius - m_innerRadius);

	if(ia > 1)
		m_aaPos[vrt] = parentCenter;
	else{
		pos_type proj, v;
		ProjectPointToRay(proj, parentCenter, m_center, m_axis);
		VecSubtract(v, parentCenter, proj);
		number len = VecLength(v);
		if(len > SMALL * avDist){	// if avDist is very small, len may be small, too
			VecScale(v, v, avDist / len);
			VecAdd(v, proj, v);

			if(ia <= 0)
				m_aaPos[vrt] = v;
			else
				VecScaleAdd(m_aaPos[vrt], 1.-ia, v, ia, parentCenter);
		}
		else
			m_aaPos[vrt] = parentCenter;
	}
}

}// end of namespace

#endif


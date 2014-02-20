// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__spherical_falloff_projector_impl__
#define __H__UG__spherical_falloff_projector_impl__

#include "spherical_falloff_projector.h"

namespace ug{

template <class TAPosition>
SphericalFalloffProjector<TAPosition>::
SphericalFalloffProjector() :
	m_pGrid(NULL)
{
}

template <class TAPosition>
SphericalFalloffProjector<TAPosition>::
SphericalFalloffProjector(Grid& grid, TAPosition& aPos,
						 const typename TAPosition::ValueType& center,
						 number innerRadius, number outerRadius) :
	m_pGrid(&grid),
	m_center(center),
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

template <class TAPosition>
SphericalFalloffProjector<TAPosition>::
~SphericalFalloffProjector()
{
}

template <class TAPosition>
void SphericalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = m_aaPos[parent];
}

template <class TAPosition>
void SphericalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	perform_projection(vrt, parent);
}

template <class TAPosition>
void SphericalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
	perform_projection(vrt, parent);
}

template <class TAPosition>
void SphericalFalloffProjector<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
	perform_projection(vrt, parent);
}

template <class TAPosition>
int SphericalFalloffProjector<TAPosition>::
current_pos(number* coordsOut, Vertex* vrt, int maxCoords)
{
	return IRefinementCallback::current_pos_helper(coordsOut, vrt, maxCoords, m_aaPos);
}

template <class TAPosition>
template <class TElem>
void SphericalFalloffProjector<TAPosition>::
perform_projection(Vertex* vrt, TElem* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
//	first calculate the average distance of corners of parent
	typename TElem::ConstVertexArray vrts = parent->vertices();
	size_t numVrts = parent->num_vertices();
	number avDist = 0;
	pos_type parentCenter;
	VecSet(parentCenter, 0);

	for(size_t i = 0; i < numVrts; ++i){
		const pos_type& p = m_aaPos[vrts[i]];
		avDist += VecDistance(p, m_center);
		VecAdd(parentCenter, parentCenter, p);
	}

	avDist /= (number)numVrts;
	VecScale(parentCenter, parentCenter, 1. / (number)numVrts);

	number ia = (avDist - m_innerRadius) / (m_outerRadius - m_innerRadius);

	if(ia > 1)
		m_aaPos[vrt] = parentCenter;
	else{
	//	calculate cylindrical projection
		pos_type cylProj;
		VecSubtract(cylProj, parentCenter, m_center);
		VecNormalize(cylProj, cylProj);
		VecScale(cylProj, cylProj, avDist);
		VecAdd(cylProj, cylProj, m_center);

		if(ia <= 0)
			m_aaPos[vrt] = cylProj;
		else
			VecScaleAdd(m_aaPos[vrt], 1.-ia, cylProj, ia, parentCenter);
	}
}

}// end of namespace

#endif

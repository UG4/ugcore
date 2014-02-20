// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__sphere_projector_impl__
#define __H__UG__sphere_projector_impl__

#include "sphere_projector.h"

namespace ug{

template <class TAPosition>
SphereProjector<TAPosition>::
SphereProjector() :
	m_pGrid(NULL)
{
}

template <class TAPosition>
SphereProjector<TAPosition>::
SphereProjector(Grid& grid, TAPosition& aPos,
						 const typename TAPosition::ValueType& center) :
	m_pGrid(&grid),
	m_center(center)
{
//	we have to make sure that aPos is attached at the grid.
//	This is important to avoid crashes later on.
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	m_aaPos.access(grid, aPos);
}

template <class TAPosition>
SphereProjector<TAPosition>::
~SphereProjector()
{
}

template <class TAPosition>
void SphereProjector<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	m_aaPos[vrt] = m_aaPos[parent];
}

template <class TAPosition>
void SphereProjector<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	perform_projection(vrt, parent);
}

template <class TAPosition>
void SphereProjector<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
	perform_projection(vrt, parent);
}

template <class TAPosition>
void SphereProjector<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
	perform_projection(vrt, parent);
}

template <class TAPosition>
int SphereProjector<TAPosition>::
current_pos(number* coordsOut, Vertex* vrt, int maxCoords)
{
	return IRefinementCallback::current_pos_helper(coordsOut, vrt, maxCoords, m_aaPos);
}

template <class TAPosition>
template <class TElem>
void SphereProjector<TAPosition>::
perform_projection(Vertex* vrt, TElem* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");

//	first calculate the average distance of corners of parent and the parents center
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

//	calculate projection
	pos_type cylProj;
	VecSubtract(cylProj, parentCenter, m_center);
	VecNormalize(cylProj, cylProj);
	VecScale(cylProj, cylProj, avDist);
	VecAdd(m_aaPos[vrt], cylProj, m_center);
}

}// end of namespace

#endif

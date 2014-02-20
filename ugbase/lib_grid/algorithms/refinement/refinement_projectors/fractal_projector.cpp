// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#include "fractal_projector.h"

namespace ug{


////////////////////////////////////////////////////////////////////////
RefinementCallbackFractal::
RefinementCallbackFractal()
{
}

////////////////////////////////////////////////////////////////////////
RefinementCallbackFractal::
RefinementCallbackFractal(Grid& grid, number scaleFac, APosition& aPos) :
	RefinementCallbackLinear<APosition>(grid, aPos),
	m_scaleFac(scaleFac)
{
}

////////////////////////////////////////////////////////////////////////
RefinementCallbackFractal::
~RefinementCallbackFractal()
{
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackFractal::
new_vertex(Vertex* vrt, Edge* parent)
{
//	set the vertex to the center by calling the parents method
	RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);

//	calculate the normal of the edge
	vector3 n;
	CalculateNormal(n, *m_pGrid, parent, m_aaPos);
//	first scale it to the same length as the edge, then scale it with the
//	given scaleFac
	number len = VecDistance(m_aaPos[parent->vertex(0)], m_aaPos[parent->vertex(1)]);
	VecScale (n, n, m_scaleFac * len);

//	offset the vertex by the calculated normal
	VecAdd(m_aaPos[vrt], m_aaPos[vrt], n);
}

////////////////////////////////////////////////////////////////////////
void RefinementCallbackFractal::
new_vertex(Vertex* vrt, Face* parent)
{
//	set the vertex to the center by calling the parents method
	RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);

//	calculate the normal of the face
	vector3 n;
	CalculateNormal(n, parent, m_aaPos);

//	we have to alter the length.
//	first scale it to the average length of the faces edges,
//	then scale it with the given scaleFac.
	number avLen = 0;

	size_t numVrts = parent->num_vertices();
	for(size_t i = 0; i < numVrts; ++i){
		avLen += VecDistance(m_aaPos[parent->vertex(i)],
							 m_aaPos[parent->vertex((i+1)% numVrts)]);
	}
	avLen /= (number)numVrts;
	VecScale (n, n, m_scaleFac * avLen);

//	offset the vertex by the calculated normal
	VecAdd(m_aaPos[vrt], m_aaPos[vrt], n);
}

}// end of namespace

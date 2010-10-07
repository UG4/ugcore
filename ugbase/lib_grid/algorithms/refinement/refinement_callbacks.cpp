// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m07 d08

#include <cassert>
#include "refinement_callbacks.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"

namespace ug
{

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
		RefinementCallbackLinear<APosition>::new_vertex(vrt, parent);
	}
}




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
new_vertex(VertexBase* vrt, EdgeBase* parent)
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
new_vertex(VertexBase* vrt, Face* parent)
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


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
RefinementCallbackSubdivisionLoop::
RefinementCallbackSubdivisionLoop()
{
}

RefinementCallbackSubdivisionLoop::
RefinementCallbackSubdivisionLoop(MultiGrid& mg,
								  APosition& aPos) :
	BaseClass(mg, aPos),
	m_pMG(&mg)
{
}

RefinementCallbackSubdivisionLoop::
~RefinementCallbackSubdivisionLoop()
{

}

void RefinementCallbackSubdivisionLoop::
new_vertex(VertexBase* vrt, VertexBase* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	if(is_crease_vertex(parent)){
	//	get the neighboured crease edges
		EdgeBase* nbrs[2];
		size_t numNbrs = 0;
		for(Grid::AssociatedEdgeIterator iter = m_pMG->associated_edges_begin(parent);
			iter != m_pMG->associated_edges_end(parent); ++iter)
		{
			if(is_crease_edge(*iter)){
				nbrs[numNbrs] = *iter;
				++numNbrs;
				if(numNbrs == 2)
					break;
			}
		}
		
		if(numNbrs == 2){
			pos_type& p0 = m_aaPos[GetConnectedVertex(nbrs[0], parent)];
			pos_type& p1 = m_aaPos[GetConnectedVertex(nbrs[1], parent)];
			VecScaleAdd(m_aaPos[vrt], 0.75, m_aaPos[parent], 0.125, p0, 0.125, p1);
		}
		else{
			BaseClass::new_vertex(vrt, parent);
		}
	}
	else{
		BaseClass::new_vertex(vrt, parent);
	}
}

void RefinementCallbackSubdivisionLoop::
new_vertex(VertexBase* vrt, EdgeBase* parent)
{
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	if(is_crease_edge(parent)){
		BaseClass::new_vertex(vrt, parent);
	}
	else{
	//	todo: apply loop-subdivision on outer elements
		BaseClass::new_vertex(vrt, parent);
	}
}
/*
void RefinementCallbackSubdivisionLoop::
new_vertex(VertexBase* vrt, Face* parent)
{

}
*/
bool RefinementCallbackSubdivisionLoop::
is_crease_vertex(VertexBase* vrt)
{
	return IsBoundaryVertex2D(*m_pMG, vrt);
}

bool RefinementCallbackSubdivisionLoop::
is_crease_edge(EdgeBase* edge)
{
	return IsBoundaryEdge2D(*m_pMG, edge);
}

}// end of namespace

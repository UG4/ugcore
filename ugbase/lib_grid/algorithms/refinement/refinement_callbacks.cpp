// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m07 d08

#include <cassert>
#include <algorithm>
#include "refinement_callbacks.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/subdivision/subdivision_rules_piecewise_loop.h"
using namespace std;

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
	init();
}

RefinementCallbackSubdivisionLoop::
RefinementCallbackSubdivisionLoop(MultiGrid& mg,
								  APosition& aPos) :
	BaseClass(mg, aPos),
	m_pMG(&mg)
{
	init();	
}

RefinementCallbackSubdivisionLoop::
~RefinementCallbackSubdivisionLoop()
{

}

void RefinementCallbackSubdivisionLoop::
init()
{
//	precalculate betas
//	the number is quite arbitrary here.
	size_t numPrecals = 16;
	m_betas.resize(numPrecals);
	for(size_t i = 0; i < numPrecals; ++i)
		m_betas[i] = calculate_beta(i);
}

void RefinementCallbackSubdivisionLoop::
new_vertex(VertexBase* vrt, VertexBase* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	
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
			vector3 w = subdiv.ref_even_crease_weights();
			
			VecScaleAdd(m_aaPos[vrt], w.x, m_aaPos[parent],
						w.y, p0, w.z, p1);
		}
		else{
			BaseClass::new_vertex(vrt, parent);
		}
	}
	else{
	//	perform loop subdivision on even vertices
	//	first get neighboured vertices
	//todo: replace this by a method
		size_t valence = 0;
		pos_type p;
		VecSet(p, 0);

		for(Grid::AssociatedEdgeIterator iter = m_pMG->associated_edges_begin(parent);
			iter != m_pMG->associated_edges_end(parent); ++iter)
		{
			VecAdd(p, p, m_aaPos[GetConnectedVertex(*iter, parent)]);
			++valence;
		}
		
		number centerWgt = subdiv.ref_even_inner_center_weight(valence);
		number nbrWgt = subdiv.ref_even_inner_nbr_weight(valence);
		
		VecScaleAdd(m_aaPos[vrt],
					centerWgt, m_aaPos[parent],
					nbrWgt, p);
/*		
		number beta = get_beta(valence);
	
		VecScaleAdd(m_aaPos[vrt], beta, p,
					1.0 - (number)valence * beta, m_aaPos[parent]); 
*/
	
	}
}

void RefinementCallbackSubdivisionLoop::
new_vertex(VertexBase* vrt, EdgeBase* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	
	assert(m_aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	
	if(is_crease_edge(parent)){
		vector2 wghts = subdiv.ref_odd_crease_weights();
		VecScaleAdd(m_aaPos[vrt], wghts.x, m_aaPos[parent->vertex(0)],
					wghts.y, m_aaPos[parent->vertex(1)]);
	}
	else{
	//	apply loop-subdivision on inner elements
	//	get the neighboured triangles
		Face* f[2];
		if(GetAssociatedFaces(f, *m_pMG, parent, 2) == 2){
			if(f[0]->num_vertices() == 3 && f[1]->num_vertices() == 3){
			//	the 4 vertices that are important for the calculation
				VertexBase* v[4];
				v[0] = parent->vertex(0); v[1] = parent->vertex(1);
				v[2] = GetConnectedVertex(parent, f[0]);
				v[3] = GetConnectedVertex(parent, f[1]);
				
				vector4 wghts;
				
			//	THIS PIECE OF CODE CAN LEAD TO PROBLEMS REGARDING THE LIMIT PROJECTION
			//	the intention of the following code is to guarantee smoothness of
			//	the grid even at irregular crease vertices. However, the limit projection
			//	is not suited for this kind of weighting.
				bool isCrease0 = is_crease_vertex(v[0]);
				bool isCrease1 = is_crease_vertex(v[1]);
			//	if exactly one of the two is a crease vertex, special
			//	weighting has to be performed
			//todo: this does not yet work for inner creases.
				if((isCrease0 && !isCrease1) || (!isCrease0 && isCrease1))
				{
				//	the crease vertex has to be in v[0]
					if(isCrease1)
						swap(v[0], v[1]);

				//	todo: replace this with a method call
				//	get the number of edges that are connected to v[0]
				//	todo: only check edges that are on the correct side of the crease.
					size_t valence = 0;
					for(Grid::AssociatedEdgeIterator iter = m_pMG->associated_edges_begin(v[0]);
						iter != m_pMG->associated_edges_end(v[0]); ++iter)
					{
						++valence;
					}
					
					wghts = subdiv.ref_odd_inner_weights(valence);
				}
				else{
					wghts = subdiv.ref_odd_inner_weights();
				}
				
				VecScaleAdd(m_aaPos[vrt],
							wghts.x, m_aaPos[v[0]], wghts.y, m_aaPos[v[1]],
							wghts.z, m_aaPos[v[2]], wghts.w, m_aaPos[v[3]]);
				
			}
			else
				BaseClass::new_vertex(vrt, parent);				
		}
		else
			BaseClass::new_vertex(vrt, parent);
	}
}

void RefinementCallbackSubdivisionLoop::
new_vertex(VertexBase* vrt, Face* parent)
{
//	this woul'd only be interesting for quad subdivision.
	BaseClass::new_vertex(vrt, parent);
}

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

number RefinementCallbackSubdivisionLoop::
get_beta(size_t valency)
{
	if(valency < m_betas.size())
		return m_betas[valency];
		
	return calculate_beta(valency);
}

number RefinementCallbackSubdivisionLoop::
calculate_beta(size_t valency)
{
	if(valency == 6)
		return 0.0625;
		
	if(valency > 0){
		const number tmp = 0.375 + 0.25 * cos((2.0*PI)/(number)valency);
		return (0.625 - tmp*tmp)/(number)valency;
	}

	return 0;
}

}// end of namespace

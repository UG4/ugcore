// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m08 d25


#ifndef __H__LIB_GRID__REFINEMENT_CALLBACKS_IMPL__
#define __H__LIB_GRID__REFINEMENT_CALLBACKS_IMPL__

#include <cassert>
#include <algorithm>
#include "refinement_callbacks.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/algorithms/subdivision/subdivision_rules_piecewise_loop.h"

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


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class TAPosition>
RefinementCallbackSubdivisionLoop<TAPosition>::
RefinementCallbackSubdivisionLoop()
{
}

template <class TAPosition>
RefinementCallbackSubdivisionLoop<TAPosition>::
RefinementCallbackSubdivisionLoop(MultiGrid& mg,
								  TAPosition& aPos) :
	BaseClass(mg, aPos),
	m_pMG(&mg)
{
}

template <class TAPosition>
RefinementCallbackSubdivisionLoop<TAPosition>::
~RefinementCallbackSubdivisionLoop()
{

}

template <class TAPosition>
void RefinementCallbackSubdivisionLoop<TAPosition>::
new_vertex(VertexBase* vrt, VertexBase* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;
	
	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
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
			pos_type& p0 = aaPos[GetConnectedVertex(nbrs[0], parent)];
			pos_type& p1 = aaPos[GetConnectedVertex(nbrs[1], parent)];
			vector3 w = subdiv.ref_even_crease_weights();
			
			VecScaleAdd(aaPos[vrt], w.x, aaPos[parent],
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
			VecAdd(p, p, aaPos[GetConnectedVertex(*iter, parent)]);
			++valence;
		}
		
		number centerWgt = subdiv.ref_even_inner_center_weight(valence);
		number nbrWgt = subdiv.ref_even_inner_nbr_weight(valence);
		
		VecScaleAdd(aaPos[vrt],
					centerWgt, aaPos[parent],
					nbrWgt, p);
/*		
		number beta = get_beta(valence);
	
		VecScaleAdd(m_aaPos[vrt], beta, p,
					1.0 - (number)valence * beta, m_aaPos[parent]); 
*/
	
	}
}

template <class TAPosition>
void RefinementCallbackSubdivisionLoop<TAPosition>::
new_vertex(VertexBase* vrt, EdgeBase* parent)
{
	using std::swap;
	
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;
	
	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	
	if(is_crease_edge(parent)){
		vector2 wghts = subdiv.ref_odd_crease_weights();
		VecScaleAdd(aaPos[vrt], wghts.x, aaPos[parent->vertex(0)],
					wghts.y, aaPos[parent->vertex(1)]);
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
				
				VecScaleAdd(aaPos[vrt],
							wghts.x, aaPos[v[0]], wghts.y, aaPos[v[1]],
							wghts.z, aaPos[v[2]], wghts.w, aaPos[v[3]]);
				
			}
			else
				BaseClass::new_vertex(vrt, parent);				
		}
		else
			BaseClass::new_vertex(vrt, parent);
	}
}

template <class TAPosition>
void RefinementCallbackSubdivisionLoop<TAPosition>::
new_vertex(VertexBase* vrt, Face* parent)
{
//	this woul'd only be interesting for quad subdivision.
	BaseClass::new_vertex(vrt, parent);
}

template <class TAPosition>
bool RefinementCallbackSubdivisionLoop<TAPosition>::
is_crease_vertex(VertexBase* vrt)
{
	return IsBoundaryVertex2D(*m_pMG, vrt);
}

template <class TAPosition>
bool RefinementCallbackSubdivisionLoop<TAPosition>::
is_crease_edge(EdgeBase* edge)
{
	return IsBoundaryEdge2D(*m_pMG, edge);
}

}// end of namespace

#endif

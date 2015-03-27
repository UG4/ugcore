// created by Sebastian Reiter
// s.b.reiter@gmail.com
// Feb 20, 2014

#ifndef __H__UG__loop_subdivision_projectors_impl__
#define __H__UG__loop_subdivision_projectors_impl__

#include "loop_subdivision_projectors.h"
#include "lib_grid/callbacks/basic_callbacks.h"

namespace ug{

template <class TAPosition>
SubdivisionLoopBoundaryProjector<TAPosition>::
SubdivisionLoopBoundaryProjector()
{
}

template <class TAPosition>
SubdivisionLoopBoundaryProjector<TAPosition>::
SubdivisionLoopBoundaryProjector(Grid& g,
								  TAPosition& aPos,
								  TAPosition& aTargetPos) :
	BaseClass(g, aPos)
{
//	we have to make sure that aTargetPos is attached at the grid.
//	This is important to avoid crashes later on.
	if(!g.has_vertex_attachment(aTargetPos))
		g.attach_to_vertices(aTargetPos);

	m_aaTargetPos.access(g, aTargetPos);
}

template <class TAPosition>
SubdivisionLoopBoundaryProjector<TAPosition>::
~SubdivisionLoopBoundaryProjector()
{

}

template <class TAPosition>
void SubdivisionLoopBoundaryProjector<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

	if(is_crease_vertex(parent)){
	//	get the neighboured crease edges
		Edge* nbrs[2];
		size_t numNbrs = 0;
		for(Grid::AssociatedEdgeIterator iter =
			BaseClass::m_pGrid->associated_edges_begin(parent);
			iter != BaseClass::m_pGrid->associated_edges_end(parent); ++iter)
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

			VecScaleAdd(m_aaTargetPos[vrt], w.x(), aaPos[parent],
						w.y(), p0, w.z(), p1);
		}
		else{
			m_aaTargetPos[vrt] = aaPos[parent];
		}
	}
	else{
		m_aaTargetPos[vrt] = aaPos[parent];
	}
}

template <class TAPosition>
void SubdivisionLoopBoundaryProjector<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	using std::swap;

	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

	if(is_crease_edge(parent)){
		vector2 wghts = subdiv.ref_odd_crease_weights();
		VecScaleAdd(m_aaTargetPos[vrt], wghts.x(), aaPos[parent->vertex(0)],
					wghts.y(), aaPos[parent->vertex(1)]);
	}
	else{
		m_aaTargetPos[vrt] = CalculateCenter(parent, aaPos);
	}
}

template <class TAPosition>
void SubdivisionLoopBoundaryProjector<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

	m_aaTargetPos[vrt] = CalculateCenter(parent, aaPos);
}

template <class TAPosition>
void SubdivisionLoopBoundaryProjector<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

	m_aaTargetPos[vrt] = CalculateCenter(parent, aaPos);
}

template <class TAPosition>
bool SubdivisionLoopBoundaryProjector<TAPosition>::
is_crease_vertex(Vertex* vrt)
{
	return !IsRegularSurfaceVertex(*BaseClass::m_pGrid, vrt);
	//return IsBoundaryVertex2D(*BaseClass::m_pGrid, vrt);
}

template <class TAPosition>
bool SubdivisionLoopBoundaryProjector<TAPosition>::
is_crease_edge(Edge* edge)
{
	return NumAssociatedFaces(*BaseClass::m_pGrid, edge) != 2;
	//return IsBoundaryEdge2D(*BaseClass::m_pGrid, edge);
}


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class TAPosition>
SubdivisionLoopProjector<TAPosition>::
SubdivisionLoopProjector() :
	m_cbIsCrease(ConsiderNone())
{
}

template <class TAPosition>
SubdivisionLoopProjector<TAPosition>::
SubdivisionLoopProjector(Grid& g,
						  TAPosition& aPos,
						  TAPosition& aTargetPos) :
	BaseClass(g, aPos),
	m_cbIsCrease(ConsiderNone())
{
//	we have to make sure that aTargetPos is attached at the grid.
//	This is important to avoid crashes later on.
	if(!g.has_vertex_attachment(aTargetPos))
		g.attach_to_vertices(aTargetPos);

	m_aaTargetPos.access(g, aTargetPos);
}

template <class TAPosition>
SubdivisionLoopProjector<TAPosition>::
~SubdivisionLoopProjector()
{

}

template <class TAPosition>
void SubdivisionLoopProjector<TAPosition>::
consider_as_crease_edge(Grid::edge_traits::callback cbIsCrease)
{
	m_cbIsCrease = cbIsCrease;
}

template <class TAPosition>
void SubdivisionLoopProjector<TAPosition>::
new_vertex(Vertex* vrt, Vertex* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

//	check whether the vertex lies inside a volume geometry. If it does,
//	perform linear refinement.
	Grid& g = *BaseClass::m_pGrid;
	bool volumesExist = g.num<Volume>() > 0;

	if(volumesExist){
		if(!IsBoundaryVertex3D(g, parent)){
			aaPos[vrt] = aaPos[parent];
			return;
		}
	}

	if(is_crease_vertex(parent)){
	//	get the neighboured crease edges
		Edge* nbrs[2];
		size_t numNbrs = 0;
		for(Grid::AssociatedEdgeIterator iter =
			g.associated_edges_begin(parent);
			iter != g.associated_edges_end(parent); ++iter)
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

			VecScaleAdd(m_aaTargetPos[vrt], w.x(), aaPos[parent],
						w.y(), p0, w.z(), p1);
		}
		else{
			m_aaTargetPos[vrt] = aaPos[vrt];
		}
	}
	else{
	//	perform loop subdivision on even vertices
	//	first get neighboured vertices
	//todo: replace this by a method
		size_t valence = 0;
		pos_type p;
		VecSet(p, 0);

		for(Grid::AssociatedEdgeIterator iter =
			g.associated_edges_begin(parent);
			iter != g.associated_edges_end(parent); ++iter)
		{
			if((!volumesExist) || IsBoundaryEdge3D(g, *iter)){
				VecAdd(p, p, aaPos[GetConnectedVertex(*iter, parent)]);
				++valence;
			}
		}

		number centerWgt = subdiv.ref_even_inner_center_weight(valence);
		number nbrWgt = subdiv.ref_even_inner_nbr_weight(valence);

		VecScaleAdd(m_aaTargetPos[vrt],
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
void SubdivisionLoopProjector<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	using namespace std;
	using std::swap;

	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

//	check whether the parent edge lies inside a volume geometry. If it does,
//	perform linear refinement.
	Grid& g = *BaseClass::m_pGrid;
	if(g.num<Volume>() > 0){
		if(!IsBoundaryEdge3D(g, parent)){
			aaPos[vrt] = CalculateCenter(parent, aaPos);
			return;
		}
	}

	if(is_crease_edge(parent)){
		vector2 wghts = subdiv.ref_odd_crease_weights();
		VecScaleAdd(m_aaTargetPos[vrt], wghts.x(), aaPos[parent->vertex(0)],
					wghts.y(), aaPos[parent->vertex(1)]);
	}
	else{
	//	apply loop-subdivision on inner elements
	//	get the neighboured triangles
		Face* f[2];
		int numAssociatedBndFaces = 0;
		if(g.num<Volume>() > 0){
			vector<Face*> faces;
			CollectAssociated(faces, g, parent);
			for(size_t i = 0; i < faces.size(); ++i){
				if(IsBoundaryFace3D(g, faces[i])){
					if(numAssociatedBndFaces < 2){
						f[numAssociatedBndFaces] = faces[i];
					}
					++numAssociatedBndFaces;
				}
			}
		}
		else{
			numAssociatedBndFaces = GetAssociatedFaces(f, g, parent, 2);
		}
		if(numAssociatedBndFaces == 2){
			if(f[0]->num_vertices() == 3 && f[1]->num_vertices() == 3){
			//	the 4 vertices that are important for the calculation
				Vertex* v[4];
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
					for(Grid::AssociatedEdgeIterator iter =
						g.associated_edges_begin(v[0]);
						iter != g.associated_edges_end(v[0]); ++iter)
					{
						++valence;
					}

					wghts = subdiv.ref_odd_inner_weights(valence);
				}
				else{
					wghts = subdiv.ref_odd_inner_weights();
				}

				VecScaleAdd(m_aaTargetPos[vrt],
							wghts.x(), aaPos[v[0]], wghts.y(), aaPos[v[1]],
							wghts.z(), aaPos[v[2]], wghts.w(), aaPos[v[3]]);

			}
			else
				m_aaTargetPos[vrt] = CalculateCenter(parent, aaPos);
		}
		else
			m_aaTargetPos[vrt] = CalculateCenter(parent, aaPos);
	}
}

template <class TAPosition>
void SubdivisionLoopProjector<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
//	this would only be interesting for quad subdivision.
	m_aaTargetPos[vrt] = CalculateCenter(parent, BaseClass::m_aaPos);
}

template <class TAPosition>
void SubdivisionLoopProjector<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
//	here a more elaborate scheme would be nice.
	m_aaTargetPos[vrt] = CalculateCenter(parent, BaseClass::m_aaPos);
}

template <class TAPosition>
bool SubdivisionLoopProjector<TAPosition>::
is_crease_vertex(Vertex* vrt)
{
	if(BaseClass::m_pGrid->template num<Volume>() > 0)
		return false;

	if(!IsRegularSurfaceVertex(*BaseClass::m_pGrid, vrt))
		return true;

	Grid::edge_traits::secure_container edges;
	BaseClass::m_pGrid->associated_elements(edges, vrt);

	for(size_t i = 0; i < edges.size(); ++i){
		if(m_cbIsCrease(edges[i]))
			return true;
	}

	return false;
}

template <class TAPosition>
bool SubdivisionLoopProjector<TAPosition>::
is_crease_edge(Edge* edge)
{
	if(BaseClass::m_pGrid->template num<Volume>() > 0)
		return false;
	return (NumAssociatedFaces(*BaseClass::m_pGrid, edge) != 2 ||
			m_cbIsCrease(edge));
}

}// end of namespace

#endif

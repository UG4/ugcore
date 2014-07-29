// created by mstepnie
// martin.stepniewski@gcsc.uni-frankfurt.de
// Juli 14, 2014

#ifndef __H__UG__smooth_tetgrid_subdivision_projectors_impl__
#define __H__UG__smooth_tetgrid_subdivision_projectors_impl__

#include "smooth_tetgrid_subdivision_projectors.h"

namespace ug{


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
template <class TAPosition>
SubdivisionVolumesProjector<TAPosition>::
SubdivisionVolumesProjector()
{
}

template <class TAPosition>
SubdivisionVolumesProjector<TAPosition>::
SubdivisionVolumesProjector(Grid& g,
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
SubdivisionVolumesProjector<TAPosition>::
~SubdivisionVolumesProjector()
{

}

template <class TAPosition>
void SubdivisionVolumesProjector<TAPosition>::
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

	VecSet(m_aaTargetPos[vrt], 0);

	if(!volumesExist){
		UG_THROW("No volumes included in grid.");
	}


	if(IsBoundaryVertex3D(g, parent))
	{
	//	perform loop subdivision on even vertices
	//	first get neighboured vertices
	//	todo: replace this by a method
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
	}
	else
	{

		size_t valence = 0;
		pos_type p;

		std::vector<Volume*> volumes;
		CollectAssociated(volumes, g, vrt);
		UG_LOG("volumes: " << volumes.size() << std::endl);

	//	Iterate over all associated volumes
		for(Grid::AssociatedVolumeIterator vIter = g.associated_volumes_begin(vrt); vIter != g.associated_volumes_end(vrt); ++vIter)
		{
			VecSet(p, 0);
			Volume* vol = *vIter;
			++valence;

		//	TETRAHEDRON CASE
			if(vol->reference_object_id() == ROID_TETRAHEDRON)
			{
			//	Iterate over all associated vertices inside tetrahedron

				/*
				 * Alternative iteration
				 *
				for(Grid::AssociatedEdgeIterator eIter = g.associated_edges_begin(vol); eIter != g.associated_edges_end(vol); ++eIter)
				{
					Edge* e = *eIter;
					if(GetConnectedVertex(e, vrt) != NULL)
						VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
				}
				*/

				for(size_t i = 0; i < vol->num_vertices(); ++i)
				{
					if(i != GetVertexIndex(vol, vrt))
					{
						VecAdd(p, p, aaPos[vol->vertex(i)]);
						UG_LOG("aaPos[vol->vertex(" << i << "): " << aaPos[vol->vertex(i)] << std::endl);
					}
				}

			//	TODO: refer to subdivision rules object
				number centerWgt = -1.0/16;
				number nbrWgt = 17.0/48;

				//VecScaleAdd(m_aaTargetPos[vrt], centerWgt, aaPos[parent], nbrWgt, p);
				VecScaleAppend(m_aaTargetPos[vrt], centerWgt, aaPos[parent], nbrWgt, p);
			}

		//	OCTAHEDRON CASE
			else if(vol->reference_object_id() == ROID_OCTAHEDRON)
			{
			//	Iterate over all vertices inside octahedron, first associate ones and last the opposing one
				for(Grid::AssociatedEdgeIterator eIter = g.associated_edges_begin(vol); eIter != g.associated_edges_end(vol); ++eIter)
				{
					Edge* e = *eIter;
					if(GetConnectedVertex(e, vrt) != NULL)
					{
						VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
					}
				}

				Vertex* oppVrt = vol->vertex(vol->get_opposing_object(vrt).second);

			//	TODO: refer to subdivision rules object
				number centerWgt = 3.0/8;
				number nbrWgt = 1.0/12;
				number oppNbrWgt = 7.0/24;

				//VecScaleAdd(m_aaTargetPos[vrt], centerWgt, aaPos[parent], nbrWgt, p, oppNbrWgt, aaPos[oppVrt]);
				VecScaleAppend(m_aaTargetPos[vrt], centerWgt, aaPos[parent], nbrWgt, p, oppNbrWgt, aaPos[oppVrt]);
			}

		//	UNSUPPORTED VOLUME ELEMENT CASE
			else
			{
				UG_THROW("Volume type not supported for subdivision volumes refinement.");
			}
		}

	//	Scale vertex position by the number of associated volume elements
		VecScale(m_aaTargetPos[vrt],  m_aaTargetPos[vrt], valence);
		UG_LOG("Even: " << m_aaTargetPos[vrt].x() << "; " << m_aaTargetPos[vrt].y() << "; " << m_aaTargetPos[vrt].z() << "; Val: " << valence << std::endl);
	}
}


template <class TAPosition>
void SubdivisionVolumesProjector<TAPosition>::
new_vertex(Vertex* vrt, Edge* parent)
{
	SubdivRules_PLoop& subdiv = SubdivRules_PLoop::inst();
	Grid::VertexAttachmentAccessor<TAPosition>& aaPos = BaseClass::m_aaPos;

	assert(aaPos.valid() && "make sure to initialise the refiner-callback correctly.");
	assert(m_aaTargetPos.valid() && "make sure to initialise the refiner-callback correctly.");

//	check whether the vertex lies inside a volume geometry. If it does,
//	perform linear refinement.
	Grid& g = *BaseClass::m_pGrid;
	bool volumesExist = g.num<Volume>() > 0;

	VecSet(m_aaTargetPos[vrt], 0);

//	Init vertex position by linear refinement
	aaPos[vrt] = CalculateCenter(parent, aaPos);

	if(!volumesExist){
		UG_THROW("No volumes included in grid.");
	}

	//if(IsBoundaryVertex3D(g, vrt))

	if(IsBoundaryEdge3D(g, parent))
	{
	//	apply loop-subdivision on inner elements
	//	get the neighboured triangles
		Face* f[2];
		int numAssociatedBndFaces = 0;
		if(g.num<Volume>() > 0){
			std::vector<Face*> faces;
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

				wghts = subdiv.ref_odd_inner_weights();

				VecScaleAdd(m_aaTargetPos[vrt],
							wghts.x(), aaPos[v[0]], wghts.y(), aaPos[v[1]],
							wghts.z(), aaPos[v[2]], wghts.w(), aaPos[v[3]]);

			}
			else
				UG_THROW("Non triangular faces included in grid.");
		}
		else
			UG_THROW("numAssociatedBndFaces != 2.");
	}
	else
	{
		size_t valence = 0;
		pos_type p;

	//	Iterate over all associated volumes
		for(Grid::AssociatedVolumeIterator vIter = g.associated_volumes_begin(vrt); vIter != g.associated_volumes_end(vrt); ++vIter)
		{
			VecSet(p, 0);
			Volume* vol = *vIter;
			++valence;

		//	TETRAHEDRON CASE
			if(vol->reference_object_id() == ROID_TETRAHEDRON)
			{
			//	Iterate over all associated vertices inside tetrahedron

				/*
				 * Alternative iteration
				 *
				for(Grid::AssociatedEdgeIterator eIter = g.associated_edges_begin(vol); eIter != g.associated_edges_end(vol); ++eIter)
				{
					Edge* e = *eIter;
					if(GetConnectedVertex(e, vrt) != NULL)
						VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
				}
				*/

				for(size_t i = 0; i < vol->num_vertices(); ++i)
				{
					if(i != GetVertexIndex(vol, vrt))
						VecAdd(p, p, aaPos[vol->vertex(i)]);
				}

			//	TODO: refer to subdivision rules object
				number centerWgt = -1.0/16;
				number nbrWgt = 17.0/48;

				//VecScaleAdd(m_aaTargetPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p);
				VecScaleAppend(m_aaTargetPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p);
			}

		//	OCTAHEDRON CASE
			else if(vol->reference_object_id() == ROID_OCTAHEDRON)
			{
			//	Iterate over all vertices inside octahedron, first associate ones and last the opposing one
				for(Grid::AssociatedEdgeIterator eIter = g.associated_edges_begin(vol); eIter != g.associated_edges_end(vol); ++eIter)
				{
					Edge* e = *eIter;
					if(GetConnectedVertex(e, vrt) != NULL)
					{
						VecAdd(p, p, aaPos[GetConnectedVertex(e, vrt)]);
					}
				}

				//Vertex oppVrt = vol->get_opposing_object(vrt);
				Vertex* oppVrt = vol->vertex(vol->get_opposing_object(vrt).second);

			//	TODO: refer to subdivision rules object
				number centerWgt = 3.0/8;
				number nbrWgt = 1.0/12;
				number oppNbrWgt = 7.0/24;

				//VecScaleAdd(m_aaTargetPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p, oppNbrWgt, aaPos[oppVrt]);
				VecScaleAppend(m_aaTargetPos[vrt], centerWgt, aaPos[vrt], nbrWgt, p, oppNbrWgt, aaPos[oppVrt]);
			}

		//	UNSUPPORTED VOLUME ELEMENT CASE
			else
			{
				UG_THROW("Volume type not supported for subdivision volumes refinement.");
			}
		}

	//	Scale vertex position by the number of associated volume elements
		VecScale(m_aaTargetPos[vrt],  m_aaTargetPos[vrt], valence);
		UG_LOG("Odd: " << m_aaTargetPos[vrt].x() << "; " << m_aaTargetPos[vrt].y() << "; " << m_aaTargetPos[vrt].z() << std::endl);
	}
}

template <class TAPosition>
void SubdivisionVolumesProjector<TAPosition>::
new_vertex(Vertex* vrt, Face* parent)
{
//	this would only be interesting for quad subdivision.
	m_aaTargetPos[vrt] = CalculateCenter(parent, BaseClass::m_aaPos);
}

template <class TAPosition>
void SubdivisionVolumesProjector<TAPosition>::
new_vertex(Vertex* vrt, Volume* parent)
{
//	here a more elaborate scheme would be nice.
	m_aaTargetPos[vrt] = CalculateCenter(parent, BaseClass::m_aaPos);
}

template <class TAPosition>
bool SubdivisionVolumesProjector<TAPosition>::
is_crease_vertex(Vertex* vrt)
{
	if(BaseClass::m_pGrid->template num<Volume>() > 0)
		return false;
	return !IsRegularSurfaceVertex(*BaseClass::m_pGrid, vrt);
	//return IsBoundaryVertex2D(*BaseClass::m_pGrid, vrt);
}

template <class TAPosition>
bool SubdivisionVolumesProjector<TAPosition>::
is_crease_edge(Edge* edge)
{
	if(BaseClass::m_pGrid->template num<Volume>() > 0)
		return false;
	return NumAssociatedFaces(*BaseClass::m_pGrid, edge) != 2;
	//return IsBoundaryEdge2D(*BaseClass::m_pGrid, edge);
}

}// end of namespace

#endif

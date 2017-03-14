/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#include <vector>
#include "hanging_node_refiner_grid.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "ref_mark_adjusters/aniso_mark_adjuster.h"

using namespace std;

namespace ug{

HangingNodeRefiner_Grid::
HangingNodeRefiner_Grid(SPRefinementProjector projector) :
	BaseClass(projector),
	m_pGrid(NULL),
	m_aVertex(false)
{

}

HangingNodeRefiner_Grid::
HangingNodeRefiner_Grid(Grid& grid,
						SPRefinementProjector projector) :
	BaseClass(projector),
	m_pGrid(NULL),
	m_aVertex(false)
{
	set_grid(&grid);
}

HangingNodeRefiner_Grid::
~HangingNodeRefiner_Grid()
{
	set_grid(NULL);
}

void HangingNodeRefiner_Grid::
grid_to_be_destroyed(Grid* grid)
{
	set_grid(NULL);
}

void HangingNodeRefiner_Grid::
assign_grid(Grid& grid)
{
	set_grid(&grid);
}

void HangingNodeRefiner_Grid::
set_grid(Grid* grid)
{
//	call base class implementation
	BaseClass::set_grid(grid);

//	clear own attachments
	if(m_pGrid)
	{
		m_pGrid->detach_from_edges(m_aVertex);
		m_pGrid->detach_from_faces(m_aVertex);
		m_pGrid->detach_from_faces(m_aAnisoMark);
		m_pGrid->detach_from_volumes(m_aAnisoMark);
		m_aaAnisoMark.invalidate();
		m_pGrid = NULL;
	}

//	attach new attachments
	if(grid){
		grid->attach_to_edges_dv(m_aVertex, NULL, false);
		grid->attach_to_faces_dv(m_aVertex, NULL, false);

		m_aaVertexEDGE.access(*grid, m_aVertex);
		m_aaVertexFACE.access(*grid, m_aVertex);

		m_pGrid = grid;
	}
}

bool HangingNodeRefiner_Grid::mark(Vertex* v, RefinementMark refMark)
{
	if(refMark != RM_COARSEN)
		return BaseClass::mark(v, refMark);
	return false;
}

bool HangingNodeRefiner_Grid::mark(Edge* e, RefinementMark refMark)
{
	if(refMark != RM_COARSEN)
		return BaseClass::mark(e, refMark);
	return false;
}

bool HangingNodeRefiner_Grid::mark(Face* f, RefinementMark refMark)
{
	if(refMark != RM_COARSEN)
		return BaseClass::mark(f, refMark);
	return false;
}

bool HangingNodeRefiner_Grid::mark(Volume* v, RefinementMark refMark)
{
	if(refMark != RM_COARSEN)
		return BaseClass::mark(v, refMark);
	return false;
}


void HangingNodeRefiner_Grid::
mark_aniso(Face* e, int anisoMark)
{
	if(!m_aaAnisoMark.is_valid_face_accessor())
		attach_aniso_marks();
	
	m_aaAnisoMark[e] = anisoMark;
	mark(e, RM_ANISOTROPIC);
}

void HangingNodeRefiner_Grid::
mark_aniso(Volume* e, int anisoMark)
{
	if(!m_aaAnisoMark.is_valid_volume_accessor())
		attach_aniso_marks();
	
	m_aaAnisoMark[e] = anisoMark;
	mark(e, RM_ANISOTROPIC);
}


int HangingNodeRefiner_Grid::
get_aniso_mark(Face* e) const
{
	if(m_aaAnisoMark.is_valid_face_accessor())
		return m_aaAnisoMark[e];
	return 0;
}

int HangingNodeRefiner_Grid::
get_aniso_mark(Volume* e) const
{
	if(m_aaAnisoMark.is_valid_volume_accessor())
		return m_aaAnisoMark[e];
	return 0;
}


void HangingNodeRefiner_Grid::
attach_aniso_marks()
{
//todo: more elaborate check that the method is not executed multiple times per object
	if(!m_aaAnisoMark.is_valid_face_accessor()){
		add_ref_mark_adjuster(AnisoMarkAdjuster::create());
		m_pGrid->attach_to_faces_dv(m_aAnisoMark, 0);
		m_pGrid->attach_to_volumes_dv(m_aAnisoMark, 0);
		m_aaAnisoMark.access(*m_pGrid, m_aAnisoMark, false, false, true, true);
	}
}


void HangingNodeRefiner_Grid::
num_marked_edges_local(std::vector<int>& numMarkedEdgesOut)
{
	num_marked_elems<Edge>(numMarkedEdgesOut);
}

void HangingNodeRefiner_Grid::
num_marked_faces_local(std::vector<int>& numMarkedFacesOut)
{
	num_marked_elems<Face>(numMarkedFacesOut);
}

void HangingNodeRefiner_Grid::
num_marked_volumes_local(std::vector<int>& numMarkedVolsOut)
{
	num_marked_elems<Volume>(numMarkedVolsOut);
}


template <class TElem>
void HangingNodeRefiner_Grid::
num_marked_elems(std::vector<int>& numMarkedElemsOut)
{
	numMarkedElemsOut.clear();
	numMarkedElemsOut.resize(1, 0);
	if(m_pGrid)
		numMarkedElemsOut[0] = get_refmark_selector().num<TElem>();
}

void HangingNodeRefiner_Grid::
collect_objects_for_refine()
{
	BaseClass::collect_objects_for_refine();

	Selector& sel = get_refmark_selector();

	for(ConstrainingEdgeIterator iter = sel.begin<ConstrainingEdge>();
		iter != sel.end<ConstrainingEdge>(); ++iter)
	{
		ConstrainingEdge* cge = *iter;
		if(cge->num_constrained_vertices() > 0)
			set_center_vertex(cge, cge->constrained_vertex(0));
	}

	for(ConstrainingTriangleIterator iter = sel.begin<ConstrainingTriangle>();
		iter != sel.end<ConstrainingTriangle>(); ++iter)
	{
		ConstrainingTriangle* cgt = *iter;
		if(cgt->num_constrained_vertices() > 0)
			set_center_vertex(cgt, cgt->constrained_vertex(0));
	}

	for(ConstrainingQuadrilateralIterator iter = sel.begin<ConstrainingQuadrilateral>();
		iter != sel.end<ConstrainingQuadrilateral>(); ++iter)
	{
		ConstrainingQuadrilateral* cgq = *iter;
		if(cgq->num_constrained_vertices() > 0)
			set_center_vertex(cgq, cgq->constrained_vertex(0));
	}
}

/* pre-refine
//	Resize the attachment containers
	{
		Selector& sel = get_refmark_selector();

		HNODE_PROFILE_BEGIN("HNode_ReserveAttachmentMemory");

		HNODE_PROFILE_BEGIN(HNODE_ReserveVrtData);
		mg.reserve<Vertex>(grid.num<Vertex>() +
					+ sel.num<Vertex>() + sel.num<Edge>()
					+ sel.num<Quadrilateral>() + sel.num<Hexahedron>());
		HNODE_PROFILE_END();

		HNODE_PROFILE_BEGIN(HNODE_ReserveEdgeData);
		mg.reserve<Edge>(mg.num<Edge>()
					+ 2 * mg.num<Edge>() + 3 * mg.num<Triangle>()
					+ 4 * mg.num<Quadrilateral>() + 3 * mg.num<Prism>()
					+ mg.num<Tetrahedron>()
					+ 4 * mg.num<Pyramid>() + 6 * mg.num<Hexahedron>());
		HNODE_PROFILE_END();

		HNODE_PROFILE_BEGIN(HNODE_ReserveFaceData);
		mg.reserve<Face>(mg.num<Face>()
					+ 4 * mg.num<Face>(l) + 10 * mg.num<Prism>(l)
					+ 8 * mg.num<Tetrahedron>(l)
					+ 9 * mg.num<Pyramid>(l) + 12 * mg.num<Hexahedron>(l));
		HNODE_PROFILE_END();

		HNODE_PROFILE_BEGIN(HNODE_ReserveVolData);
		mg.reserve<Volume>(mg.num<Volume>()
					+ 8 * mg.num<Tetrahedron>(l) + 8 * mg.num<Prism>(l)
					+ 6 * mg.num<Pyramid>(l) + 8 * mg.num<Hexahedron>(l));
		HNODE_PROFILE_END();

		HNODE_PROFILE_END();
	}
 */
void HangingNodeRefiner_Grid::
post_refine()
{
	if(!m_pGrid)
		throw(UGError("HangingNodeRefiner_Grid::post_refine: No grid assigned."));

//	erase unused elements
	UG_DLOG(LIB_GRID, 1, "  erasing elements.\n");

	Grid& grid = *m_pGrid;
	vector<Face*>	 	vFaces;
	vector<Volume*>		vVols;

//	erase faces that are no longer needed.
	if(grid.num_volumes() > 0)
	{
		FaceIterator iter = m_selMarkedElements.begin<Face>();
		while(iter != m_selMarkedElements.end<Face>())
		{
			Face* f = *iter;
			++iter;
			CollectVolumes(vVols, grid, f);
			if(vVols.empty())
			{
			//	erase
				grid.erase(f);
			}
		}
	}

//	erase edges that are no longer needed.
	if(grid.num_faces() > 0)
	{
		EdgeIterator iter = m_selMarkedElements.begin<Edge>();
		while(iter != m_selMarkedElements.end<Edge>())
		{
			Edge* e = *iter;
			++iter;
			CollectFaces(vFaces, grid, e);
			if(vFaces.empty())
			{
			//	erase
				grid.erase(e);
			}
		}
	}
}

void HangingNodeRefiner_Grid::
process_constraining_edge(ConstrainingEdge* e)
{
//	call original implementation
	BaseClass::process_constraining_edge(e);

//	if there are no faces, the edge can be erased
	if(m_pGrid->num_faces() == 0)
		m_pGrid->erase(e);
}

void HangingNodeRefiner_Grid::
refine_edge_with_normal_vertex(Edge* e, Vertex** newCornerVrts)
{
//	call original implementation
	BaseClass::refine_edge_with_normal_vertex(e, newCornerVrts);

//	if there are no faces, the edge can be erased
	if(m_pGrid->num_faces() == 0)
		m_pGrid->erase(e);
}

void HangingNodeRefiner_Grid::
refine_face_with_normal_vertex(Face* f, Vertex** newCornerVrts)
{
//	call original implementation
	BaseClass::refine_face_with_normal_vertex(f, newCornerVrts);

//	if there are no volumes, the face can be erased
	if(m_pGrid->num_volumes() == 0)
		m_pGrid->erase(f);
}

void HangingNodeRefiner_Grid::
process_constraining_face(ConstrainingFace* f)
{
//	call original implementation
	BaseClass::process_constraining_face(f);

//	if there are no volumes, the face can be erased
	if(m_pGrid->num_volumes() == 0)
		m_pGrid->erase(f);
}

void HangingNodeRefiner_Grid::
refine_volume_with_normal_vertex(Volume* v, Vertex** newCornerVrts)
{
//	call original implementation
	BaseClass::refine_volume_with_normal_vertex(v, newCornerVrts);

//	erase the volume
	m_pGrid->erase(v);
}

Vertex* HangingNodeRefiner_Grid::
get_center_vertex(Edge* e)
{
	return m_aaVertexEDGE[e];
}


void HangingNodeRefiner_Grid::
set_center_vertex(Edge* e, Vertex* v)
{
	m_aaVertexEDGE[e] = v;
}


Vertex* HangingNodeRefiner_Grid::
get_center_vertex(Face* f)
{
	return m_aaVertexFACE[f];
}


void HangingNodeRefiner_Grid::
set_center_vertex(Face* f, Vertex* v)
{
	m_aaVertexFACE[f] = v;
}

}// end of namespace

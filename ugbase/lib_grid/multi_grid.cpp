/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

#include <algorithm>
#include "multi_grid.h"
#include "lib_grid_messages.h"

using namespace std;

namespace ug
{

MultiGrid::MultiGrid() :
	Grid(),
	m_aVertexInfo("MultiGrid_VertexInfo"),
	m_aEdgeInfo("MultiGrid_EdgeInfo"),
	m_aFaceInfo("MultiGrid_FaceInfo"),
	m_aVolumeInfo("MultiGrid_VolumeInfo"),
	m_aParentType("MultiGrid_ParentType")
{
	init();
}

MultiGrid::MultiGrid(uint options) :
	Grid(options),
	m_aVertexInfo("MultiGrid_VertexInfo"),
	m_aEdgeInfo("MultiGrid_EdgeInfo"),
	m_aFaceInfo("MultiGrid_FaceInfo"),
	m_aVolumeInfo("MultiGrid_VolumeInfo"),
	m_aParentType("MultiGrid_ParentType")
{
	init();
}

MultiGrid::~MultiGrid()
{
//	release child infos
	for(FaceIterator iter = begin<Face>(); iter != end<Face>(); ++iter)
		release_child_info(*iter);

	for(VolumeIterator iter = begin<Volume>(); iter != end<Volume>(); ++iter)
		release_child_info(*iter);
}

void MultiGrid::init()
{
//	the subset-handler that manages the hierarchy
//	has to be registered before the multi-grid (order of create-methods).
	m_hierarchy.assign_grid(*this);
	m_hierarchy.enable_subset_inheritance(false);
	m_bHierarchicalInsertion = true;

//	the MultiGrid observes itself (its underlying grid).
	register_observer(this, OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
							OT_FACE_OBSERVER | OT_VOLUME_OBSERVER);

//	attach parent-pointers
	attach_to_faces(m_aParent);
	attach_to_volumes(m_aParent);

//	attach elem-infos
	attach_to_vertices(m_aVertexInfo);
	attach_to_edges(m_aEdgeInfo);
	attach_to_faces_dv(m_aFaceInfo, nullptr);
	attach_to_volumes_dv(m_aVolumeInfo, nullptr);

	attach_to_all(m_aParentType);

//	init accessors
	m_aaVrtInf.access(*this, m_aVertexInfo);
	m_aaEdgeInf.access(*this, m_aEdgeInfo);
	m_aaFaceInf.access(*this, m_aFaceInfo);
	m_aaParentFACE.access(*this, m_aParent);
	m_aaVolInf.access(*this, m_aVolumeInfo);
	m_aaParentVOL.access(*this, m_aParent);

	m_aaParentType.access(*this, m_aParentType, true, true, true, true);
}

void MultiGrid::create_levels(int numLevels)
{
	for(int i = 0; i < numLevels; ++i){
	//	inform the hierarchy handler, that one level has to be added
		m_hierarchy.subset_required(num_levels());
	//	send a message, that a new level has been created
		message_hub()->post_message(
				GridMessage_MultiGridChanged(GMMGCT_LEVEL_ADDED, num_levels()));
	}
}

void MultiGrid::enable_hierarchical_insertion(bool bEnable)
{
	m_bHierarchicalInsertion = bEnable;
}

////////////////////////////////////////////////////////////////////////
//	create methods
VertexIterator MultiGrid::
create_by_cloning(Vertex* pCloneMe, int level)
{
	VertexIterator iter = Grid::create_by_cloning(pCloneMe);
//	put the element into the hierarchy
//	(by default it already was assigned to level 0)
	if(level > 0){
		level_required(level);
		m_hierarchy.assign_subset(*iter, level);
	}
	return iter;
}

EdgeIterator MultiGrid::
create_by_cloning(Edge* pCloneMe, const EdgeVertices& ev, int level)
{
	EdgeIterator iter = Grid::create_by_cloning(pCloneMe, ev);
//	put the element into the hierarchy
//	(by default it already was assigned to level 0)
	if(level > 0){
		level_required(level);
		m_hierarchy.assign_subset(*iter, level);
	}
	return iter;
}

FaceIterator MultiGrid::
create_by_cloning(Face* pCloneMe, const FaceVertices& fv, int level)
{
	FaceIterator iter = Grid::create_by_cloning(pCloneMe, fv);
//	put the element into the hierarchy
//	(by default it already was assigned to level 0)
	if(level > 0){
		level_required(level);
		m_hierarchy.assign_subset(*iter, level);
	}
	return iter;
}

VolumeIterator MultiGrid::
create_by_cloning(Volume* pCloneMe, const VolumeVertices& vv, int level)
{
	VolumeIterator iter = Grid::create_by_cloning(pCloneMe, vv);
//	put the element into the hierarchy
//	(by default it already was assigned to level 0)
	if(level > 0){
		level_required(level);
		m_hierarchy.assign_subset(*iter, level);
	}
	return iter;
}


GridObject* MultiGrid::get_parent(GridObject* parent) const
{
	int baseType = parent->base_object_id();
	switch(baseType)
	{
		case VERTEX:	return get_parent((Vertex*)parent);
		case EDGE:		return get_parent((Edge*)parent);
		case FACE:		return get_parent((Face*)parent);
		case VOLUME:	return get_parent((Volume*)parent);
	}
	return nullptr;
}

////////////////////////////////////////////////////////////////////////
//	grid-observer callbacks

void MultiGrid::elements_to_be_cleared(Grid* grid)
{
//TODO: runtime of clear can be optimized in this method.
}

//	vertices
void MultiGrid::vertex_created(Grid* grid, Vertex* vrt,
								GridObject* pParent,
								bool replacesParent)
{
//	if hierarchical_insertion is disabled, the elemenet is inserted
//	into the same level as its parent.
//	From the standpoint of a multigrid-hierarchy it thus makes sense
//	to make pParents parent the parent of elem!!!

	if(replacesParent){
	//	the object given in parent will be replaced be the newly created one.
	//	The parent of pParent is thus the real parent of the new object.
		UG_ASSERT(pParent, "A parent has to exist if it shall be replaced.");
		UG_ASSERT(pParent->base_object_id() == VERTEX,
				  "only objects of the same base type can be replaced.");
		Vertex* pReplaceMe = static_cast<Vertex*>(pParent);
		GridObject* realParent = get_parent(pReplaceMe);

	//	we call a version of element_created, which allows a replace
		if(realParent){
			int baseType = realParent->base_object_id();
			switch(baseType)
			{
			case VERTEX:	element_created(vrt, (Vertex*)realParent, pReplaceMe); break;
			case EDGE:		element_created(vrt, (Edge*)realParent, pReplaceMe); break;
			case FACE:		element_created(vrt, (Face*)realParent, pReplaceMe); break;
			case VOLUME:	element_created(vrt, (Volume*)realParent, pReplaceMe); break;
			}
		}
		else
			element_created<Vertex, Vertex>(vrt, nullptr, pReplaceMe);

	//	copy pReplaceMe's children and replace parent of children
		MGVertexInfo& myInfo = get_info(vrt);
		MGVertexInfo& replaceInfo = get_info(pReplaceMe);

		if(replaceInfo.child_vertex()){
			myInfo.add_child(replaceInfo.child_vertex());
			set_parent(replaceInfo.child_vertex(), vrt);
			set_parent_type(replaceInfo.child_vertex(), VERTEX);
		}
	}
	else{
		if(!hierarchical_insertion_enabled() && pParent)
			pParent = get_parent(pParent);

		if(pParent)
		{
			int baseType = pParent->base_object_id();
			switch(baseType)
			{
			case VERTEX:	element_created(vrt, (Vertex*)pParent); break;
			case EDGE:		element_created(vrt, (Edge*)pParent); break;
			case FACE:		element_created(vrt, (Face*)pParent); break;
			case VOLUME:	element_created(vrt, (Volume*)pParent); break;
			}
		}
		else
			element_created(vrt);
	}
}

void MultiGrid::vertex_to_be_erased(Grid* grid, Vertex* vrt,
									 Vertex* replacedBy)
{
//	if replacedBy != nullptr, then vertex_created already handled the
//	deregistration at the parent.
	if(replacedBy)
		return;

	GridObject* pParent = get_parent(vrt);
	if(pParent)
	{
		int baseType = pParent->base_object_id();
		switch(baseType)
		{
		case VERTEX:	element_to_be_erased(vrt, (Vertex*)pParent); break;
		case EDGE:		element_to_be_erased(vrt, (Edge*)pParent); break;
		case FACE:		element_to_be_erased(vrt, (Face*)pParent); break;
		case VOLUME:	element_to_be_erased(vrt, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(vrt);
}

//	edges
void MultiGrid::edge_created(Grid* grid, Edge* edge,
							GridObject* pParent,
							bool replacesParent)
{
	if(replacesParent){
	//	the object given in parent will be replaced be the newly created one.
	//	The parent of pParent is thus the real parent of the new object.
		UG_ASSERT(pParent, "A parent has to exist if it shall be replaced.");
		UG_ASSERT(pParent->base_object_id() == EDGE,
				  "only objects of the same base type can be replaced.");
		Edge* pReplaceMe = static_cast<Edge*>(pParent);
		GridObject* realParent = get_parent(pReplaceMe);
		if(realParent){
		//	we call a version of element_created, which allows a replace
			int baseType = realParent->base_object_id();
			switch(baseType)
			{
			case EDGE:		element_created(edge, (Edge*)realParent, pReplaceMe); break;
			case FACE:		element_created(edge, (Face*)realParent, pReplaceMe); break;
			case VOLUME:	element_created(edge, (Volume*)realParent, pReplaceMe); break;
			}
		}
		else
			element_created<Edge, Edge>(edge, nullptr, pReplaceMe);

	//	copy pReplaceMes children and replace parent of children
		MGEdgeInfo& myInfo = get_info(edge);
		MGEdgeInfo& replaceInfo = get_info(pReplaceMe);

		if(replaceInfo.child_vertex()){
			myInfo.add_child(replaceInfo.child_vertex());
			set_parent(replaceInfo.child_vertex(), edge);
			set_parent_type(replaceInfo.child_vertex(), EDGE);
		}

		for(size_t i = 0; i < replaceInfo.num_child_edges(); ++i){
			myInfo.add_child(replaceInfo.child_edge(i));
			set_parent(replaceInfo.child_edge(i), edge);
			set_parent_type(replaceInfo.child_edge(i), EDGE);
		}
	}
	else{
		if(!hierarchical_insertion_enabled() && pParent)
			pParent = get_parent(pParent);

		if(pParent)
		{
			int baseType = pParent->base_object_id();
			switch(baseType)
			{
			case EDGE:		element_created(edge, (Edge*)pParent); break;
			case FACE:		element_created(edge, (Face*)pParent); break;
			case VOLUME:	element_created(edge, (Volume*)pParent); break;
			}
		}
		else
			element_created(edge);
	}
}

void MultiGrid::edge_to_be_erased(Grid* grid, Edge* edge,
									Edge* replacedBy)
{
	if(replacedBy)
		return;

	GridObject* pParent = get_parent(edge);
	if(pParent)
	{
		int baseType = pParent->base_object_id();
		switch(baseType)
		{
		case EDGE:		element_to_be_erased(edge, (Edge*)pParent); break;
		case FACE:		element_to_be_erased(edge, (Face*)pParent); break;
		case VOLUME:	element_to_be_erased(edge, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(edge);
}

//	faces
void MultiGrid::face_created(Grid* grid, Face* face,
							GridObject* pParent,
							bool replacesParent)
{
	if(replacesParent){
	//	the object given in parent will be replaced be the newly created one.
	//	The parent of pParent is thus the real parent of the new object.
		UG_ASSERT(pParent, "A parent has to exist if it shall be replaced.");
		UG_ASSERT(pParent->base_object_id() == FACE,
				  "only objects of the same base type can be replaced.");
		Face* pReplaceMe = static_cast<Face*>(pParent);
		GridObject* realParent = get_parent(pReplaceMe);

	//	we call a version of element_created, which allows a replace
		if(realParent){
			int baseType = realParent->base_object_id();
			switch(baseType)
			{
			case FACE:		element_created(face, (Face*)realParent, pReplaceMe); break;
			case VOLUME:	element_created(face, (Volume*)realParent, pReplaceMe); break;
			}
		}
		else
			element_created<Face, Face>(face, nullptr, pReplaceMe);

	//	copy pReplaceMes children and replace parent of children
		if(has_children(pReplaceMe)){
			create_child_info(face);
			MGFaceInfo& myInfo = get_info(face);
			MGFaceInfo& replaceInfo = get_info(pReplaceMe);

			if(replaceInfo.child_vertex()){
				myInfo.add_child(replaceInfo.child_vertex());
				set_parent(replaceInfo.child_vertex(), face);
				set_parent_type(replaceInfo.child_vertex(), FACE);
			}

			for(size_t i = 0; i < replaceInfo.num_child_edges(); ++i){
				myInfo.add_child(replaceInfo.child_edge(i));
				set_parent(replaceInfo.child_edge(i), face);
				set_parent_type(replaceInfo.child_edge(i), FACE);
			}

			for(size_t i = 0; i < replaceInfo.num_child_faces(); ++i){
				myInfo.add_child(replaceInfo.child_face(i));
				set_parent(replaceInfo.child_face(i), face);
				set_parent_type(replaceInfo.child_face(i), FACE);
			}
		}
	}
	else{
		if(!hierarchical_insertion_enabled() && pParent)
			pParent = get_parent(pParent);

		if(pParent)
		{
			int baseType = pParent->base_object_id();
			switch(baseType)
			{
			case FACE:		element_created(face, (Face*)pParent); break;
			case VOLUME:	element_created(face, (Volume*)pParent); break;
			}
		}
		else
			element_created(face);
	}
}

void MultiGrid::face_to_be_erased(Grid* grid, Face* face,
								 Face* replacedBy)
{
	if(replacedBy)
		return;

	GridObject* pParent = get_parent(face);
	if(pParent)
	{
		int baseType = pParent->base_object_id();
		switch(baseType)
		{
		case FACE:		element_to_be_erased(face, (Face*)pParent); break;
		case VOLUME:	element_to_be_erased(face, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(face);
}

//	volumes
void MultiGrid::volume_created(Grid* grid, Volume* vol,
								GridObject* pParent,
								bool replacesParent)
{
	if(replacesParent){
	//	the object given in parent will be replaced be the newly created one.
	//	The parent of pParent is thus the real parent of the new object.
		UG_ASSERT(pParent, "A parent has to exist if it shall be replaced.");
		UG_ASSERT(pParent->base_object_id() == VOLUME,
				  "only objects of the same base type can be replaced.");
		Volume* pReplaceMe = static_cast<Volume*>(pParent);
		GridObject* realParent = get_parent(pReplaceMe);

	//	we call a version of element_created, which allows a replace
		element_created(vol, (Volume*)realParent, pReplaceMe);

	//	copy pReplaceMes children and replace parent of children
		if(has_children(pReplaceMe)){
			create_child_info(vol);
			MGVolumeInfo& myInfo = get_info(vol);
			MGVolumeInfo& replaceInfo = get_info(pReplaceMe);

			if(replaceInfo.child_vertex()){
				myInfo.add_child(replaceInfo.child_vertex());
				set_parent(replaceInfo.child_vertex(), vol);
				set_parent_type(replaceInfo.child_vertex(), VOLUME);
			}

			for(size_t i = 0; i < replaceInfo.num_child_edges(); ++i){
				myInfo.add_child(replaceInfo.child_edge(i));
				set_parent(replaceInfo.child_edge(i), vol);
				set_parent_type(replaceInfo.child_edge(i), VOLUME);
			}

			for(size_t i = 0; i < replaceInfo.num_child_faces(); ++i){
				myInfo.add_child(replaceInfo.child_face(i));
				set_parent(replaceInfo.child_face(i), vol);
				set_parent_type(replaceInfo.child_face(i), VOLUME);
			}

			for(size_t i = 0; i < replaceInfo.num_child_volumes(); ++i){
				myInfo.add_child(replaceInfo.child_volume(i));
				set_parent(replaceInfo.child_volume(i), vol);
				set_parent_type(replaceInfo.child_volume(i), VOLUME);
			}
		}
	}
	else{
		if(!hierarchical_insertion_enabled() && pParent)
			pParent = get_parent(pParent);

		if(pParent)
		{
			UG_ASSERT(pParent->base_object_id() == VOLUME,
				  "Only volumes can be parents to volumes.");
			element_created(vol, (Volume*)pParent);
		}
		else
			element_created(vol);
	}
}

void MultiGrid::volume_to_be_erased(Grid* grid, Volume* vol,
									 Volume* replacedBy)
{
	if(replacedBy)
		return;

	GridObject* pParent = get_parent(vol);
	if(pParent)
	{
		UG_ASSERT(pParent->base_object_id() == VOLUME,
				  "Only volumes can be parents to volumes.");
		element_to_be_erased(vol, (Volume*)pParent);
	}
	else
		element_to_be_erased(vol);
}


void MultiGrid::check_edge_elem_infos(int level) const
{
//	check the max fill rates of each child list.
	size_t maxChildEdges = 0;

	for(ConstEdgeIterator iter = begin<Edge>(level);
		iter != end<Edge>(level); ++iter)
		maxChildEdges = max(get_info(*iter).num_child_edges(), maxChildEdges);

	UG_LOG("MultiGrid: max edge child edges on level " << level << ": " << (int)maxChildEdges << endl);
}

void MultiGrid::check_face_elem_infos(int level) const
{
//	check the max fill rates of each child list.
	size_t maxChildEdges = 0;
	size_t maxChildFaces = 0;

	for(ConstFaceIterator iter = begin<Face>(level);
		iter != end<Face>(level); ++iter)
	{
		maxChildEdges = max(get_info(*iter).num_child_edges(), maxChildEdges);
		maxChildFaces = max(get_info(*iter).num_child_faces(), maxChildFaces);
	}

	UG_LOG("MultiGrid: max face child edges on level " << level << ": " << (int)maxChildEdges << endl);
	UG_LOG("MultiGrid: max face child faces on level " << level << ": " << (int)maxChildFaces << endl);
}

void MultiGrid::check_volume_elem_infos(int level) const
{
//	check the max fill rates of each child list.
	size_t maxChildEdges = 0;
	size_t maxChildFaces = 0;
	size_t maxChildVolumes = 0;

	for(ConstVolumeIterator iter = begin<Volume>(level);
		iter != end<Volume>(level); ++iter)
	{
		maxChildEdges = max(get_info(*iter).num_child_edges(), maxChildEdges);
		maxChildFaces = max(get_info(*iter).num_child_faces(), maxChildFaces);
		maxChildVolumes = max(get_info(*iter).num_child_volumes(), maxChildVolumes);
	}

	UG_LOG("MultiGrid: max volume child edges on level " << level << ": " << (int)maxChildEdges << endl);
	UG_LOG("MultiGrid: max volume child faces on level " << level << ": " << (int)maxChildFaces << endl);
	UG_LOG("MultiGrid: max volume child volumes on level " << level << ": " << (int)maxChildVolumes << endl);
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of Info-Classes
void MGVertexInfo::unregister_from_children(MultiGrid& mg)
{
	if(m_pVrtChild)
		mg.set_parent(m_pVrtChild, nullptr);
	clear_children();
}

void MGEdgeInfo::unregister_from_children(MultiGrid& mg)
{
	if(m_pVrtChild)
		mg.set_parent(m_pVrtChild, nullptr);
	for(int i = 0; i < m_numEdgeChildren; ++i)
		mg.set_parent(m_pEdgeChild[i], nullptr);
	clear_children();
}

void MGFaceInfo::unregister_from_children(MultiGrid& mg)
{
	if(m_pVrtChild)
		mg.set_parent(m_pVrtChild, nullptr);
	for(int i = 0; i < m_numEdgeChildren; ++i)
		mg.set_parent(m_pEdgeChild[i], nullptr);
	for(int i = 0; i < m_numFaceChildren; ++i)
		mg.set_parent(m_pFaceChild[i], nullptr);
	clear_children();
}

void MGVolumeInfo::unregister_from_children(MultiGrid& mg)
{
	if(m_pVrtChild)
		mg.set_parent(m_pVrtChild, nullptr);
	for(size_t i = 0; i < m_edgeChildren.size(); ++i)
		mg.set_parent(m_edgeChildren[i], nullptr);
	for(size_t i = 0; i < m_faceChildren.size(); ++i)
		mg.set_parent(m_faceChildren[i], nullptr);
	for(size_t i = 0; i < m_volumeChildren.size(); ++i)
		mg.set_parent(m_volumeChildren[i], nullptr);
	clear_children();
}

}//	end of namespace

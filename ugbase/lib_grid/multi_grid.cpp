//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d10

#include <algorithm>
#include "multi_grid.h"

using namespace std;

namespace ug
{

MultiGrid::MultiGrid()
{
	init();
}

MultiGrid::MultiGrid(uint options) : Grid(options)
{
	init();
}

MultiGrid::~MultiGrid()
{
	unregister_observer(this);
}

void MultiGrid::init()
{
//	the subset-handler that manages the hierarchy
//	has to be registered before the multi-grid (order of create-methods).
	m_hierarchy.assign_grid(*this);
	m_hierarchy.enable_subset_inheritance(false);
	m_bHierarchicalInsertion = true;

//	the MultiGrid observes itself (its underlying grid).
	register_observer(this, OT_FULL_OBSERVER);

//	attach elem-infos
	attach_to_vertices(m_aVertexInfo);
	attach_to_edges(m_aEdgeInfo);
	attach_to_faces(m_aFaceInfo);
	attach_to_volumes(m_aVolumeInfo);

//	init accessors
	m_aaVrtInf.access(*this, m_aVertexInfo);
	m_aaEdgeInf.access(*this, m_aEdgeInfo);
	m_aaFaceInf.access(*this, m_aFaceInfo);
	m_aaVolInf.access(*this, m_aVolumeInfo);
}

void MultiGrid::enable_hierarchical_insertion(bool bEnable)
{
	m_bHierarchicalInsertion = bEnable;
}

GeometricObject* MultiGrid::get_parent(GeometricObject* parent)
{
	int baseType = parent->base_object_type_id();
	switch(baseType)
	{
	case VERTEX:	return get_parent((VertexBase*)parent);
	case EDGE:		return get_parent((EdgeBase*)parent);
	case FACE:		return get_parent((Face*)parent);
	case VOLUME:	return get_parent((Volume*)parent);
	}
	return NULL;
}

void MultiGrid::set_status(VertexBase* vrt)
{
}

void MultiGrid::set_status(EdgeBase* edge)
{
}

void MultiGrid::set_status(Face* face)
{
}

void MultiGrid::set_status(Volume* vol)
{
}
	
////////////////////////////////////////////////////////////////////////
//	grid-observer callbacks

void MultiGrid::elements_to_be_cleared(Grid* grid)
{
//TODO: runtime of clear can be optimized in this method.
}

//	vertices
void MultiGrid::vertex_created(Grid* grid, VertexBase* vrt, GeometricObject* pParent)
{
//	if hierarchical_insertion is disabled, the elemenet is inserted
//	into the same level as its parent.
//	From the standpoint of a multigrid-hierarchy it thus makes sense
//	to make pParents parent the parent of elem!!!

	if(!hierarchical_insertion_enabled() && pParent)
		pParent = get_parent(pParent);

	if(pParent)
	{
		int baseType = pParent->base_object_type_id();

		switch(baseType)
		{
		case VERTEX:	element_created(vrt, (VertexBase*)pParent); break;
		case EDGE:		element_created(vrt, (EdgeBase*)pParent); break;
		case FACE:		element_created(vrt, (Face*)pParent); break;
		case VOLUME:	element_created(vrt, (Volume*)pParent); break;
		}
	}
	else
		element_created(vrt);
}

void MultiGrid::vertex_to_be_erased(Grid* grid, VertexBase* vrt)
{
	VertexInfo& info = get_info(vrt);
	GeometricObject* pParent = info.m_pParent;
	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case VERTEX:	element_to_be_erased(vrt, (VertexBase*)pParent); break;
		case EDGE:		element_to_be_erased(vrt, (EdgeBase*)pParent); break;
		case FACE:		element_to_be_erased(vrt, (Face*)pParent); break;
		case VOLUME:	element_to_be_erased(vrt, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(vrt);
}

void MultiGrid::vertex_to_be_replaced(Grid* grid, VertexBase* vrtOld, VertexBase* vrtNew)
{
//	child and parent connections are handled by create and erase.
}

//	edges
void MultiGrid::edge_created(Grid* grid, EdgeBase* edge, GeometricObject* pParent)
{
	if(!hierarchical_insertion_enabled() && pParent)
		pParent = get_parent(pParent);

	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case EDGE:		element_created(edge, (EdgeBase*)pParent); break;
		case FACE:		element_created(edge, (Face*)pParent); break;
		case VOLUME:	element_created(edge, (Volume*)pParent); break;
		}
	}
	else
		element_created(edge);
}

void MultiGrid::edge_to_be_erased(Grid* grid, EdgeBase* edge)
{
	EdgeInfo& info = get_info(edge);
	GeometricObject* pParent = info.m_pParent;
	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case EDGE:		element_to_be_erased(edge, (EdgeBase*)pParent); break;
		case FACE:		element_to_be_erased(edge, (Face*)pParent); break;
		case VOLUME:	element_to_be_erased(edge, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(edge);
}

void MultiGrid::edge_to_be_replaced(Grid* grid, EdgeBase* edgeOld, EdgeBase* edgeNew)
{
//	child and parent connections are handled by create and erase.
}

//	faces
void MultiGrid::face_created(Grid* grid, Face* face, GeometricObject* pParent)
{
	if(!hierarchical_insertion_enabled() && pParent)
		pParent = get_parent(pParent);

	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case FACE:		element_created(face, (Face*)pParent); break;
		case VOLUME:	element_created(face, (Volume*)pParent); break;
		}
	}
	else
		element_created(face);
}

void MultiGrid::face_to_be_erased(Grid* grid, Face* face)
{
	FaceInfo& info = get_info(face);
	GeometricObject* pParent = info.m_pParent;
	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case FACE:		element_to_be_erased(face, (Face*)pParent); break;
		case VOLUME:	element_to_be_erased(face, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(face);
}

void MultiGrid::face_to_be_replaced(Grid* grid, Face* faceOld, Face* faceNew)
{
//	child and parent connections are handled by create and erase.
}

//	volumes
void MultiGrid::volume_created(Grid* grid, Volume* vol, GeometricObject* pParent)
{
	if(!hierarchical_insertion_enabled() && pParent)
		pParent = get_parent(pParent);

	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case VOLUME:	element_created(vol, (Volume*)pParent); break;
		}
	}
	else
		element_created(vol);
}

void MultiGrid::volume_to_be_erased(Grid* grid, Volume* vol)
{
	VolumeInfo& info = get_info(vol);
	GeometricObject* pParent = info.m_pParent;
	if(pParent)
	{
		int baseType = pParent->base_object_type_id();
		switch(baseType)
		{
		case VOLUME:	element_to_be_erased(vol, (Volume*)pParent); break;
		}
	}
	else
		element_to_be_erased(vol);
}

void MultiGrid::volume_to_be_replaced(Grid* grid, Volume* volOld, Volume* volNew)
{
//	child and parent connections are handled by create and erase.
}

void MultiGrid::check_edge_elem_infos(int level)
{
//	check the max fill rates of each child list.
	byte maxChildEdges = 0;

	for(EdgeBaseIterator iter = begin<EdgeBase>(level);
		iter != end<EdgeBase>(level); ++iter)
		maxChildEdges = max(get_info(*iter).m_numEdgeChildren, maxChildEdges);

	cout << "MultiGrid: max edge child edges on level " << level << ": " << (int)maxChildEdges << endl;
}

void MultiGrid::check_face_elem_infos(int level)
{
//	check the max fill rates of each child list.
	byte maxChildEdges = 0;
	byte maxChildFaces = 0;

	for(FaceIterator iter = begin<Face>(level);
		iter != end<Face>(level); ++iter)
	{
		maxChildEdges = max(get_info(*iter).m_numEdgeChildren, maxChildEdges);
		maxChildFaces = max(get_info(*iter).m_numFaceChildren, maxChildFaces);
	}

	cout << "MultiGrid: max face child edges on level " << level << ": " << (int)maxChildEdges << endl;
	cout << "MultiGrid: max face child faces on level " << level << ": " << (int)maxChildFaces << endl;
}

void MultiGrid::check_volume_elem_infos(int level)
{
//	check the max fill rates of each child list.
	byte maxChildEdges = 0;
	byte maxChildFaces = 0;
	byte maxChildVolumes = 0;

	for(VolumeIterator iter = begin<Volume>(level);
		iter != end<Volume>(level); ++iter)
	{
		maxChildEdges = max(get_info(*iter).m_numEdgeChildren, maxChildEdges);
		maxChildFaces = max(get_info(*iter).m_numFaceChildren, maxChildFaces);
		maxChildVolumes = max(get_info(*iter).m_numVolChildren, maxChildVolumes);
	}

	cout << "MultiGrid: max volume child edges on level " << level << ": " << (int)maxChildEdges << endl;
	cout << "MultiGrid: max volume child faces on level " << level << ": " << (int)maxChildFaces << endl;
	cout << "MultiGrid: max volume child volumes on level " << level << ": " << (int)maxChildVolumes << endl;
}



////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	implementation of Info-Classes
void MGVertexInfo::erase_all_children(Grid& grid)
{
	if(m_pVrtChild)
		grid.erase(m_pVrtChild);
	clear();
}

void MGEdgeInfo::erase_all_children(Grid& grid)
{
	if(m_pVrtChild)
		grid.erase(m_pVrtChild);
	for(int i = 0; i < m_numEdgeChildren; ++i)
		grid.erase(m_pEdgeChild[i]);
	clear();
}

void MGFaceInfo::erase_all_children(Grid& grid)
{
	if(m_pVrtChild)
		grid.erase(m_pVrtChild);
	for(int i = 0; i < m_numEdgeChildren; ++i)
		grid.erase(m_pEdgeChild[i]);
	for(int i = 0; i < m_numFaceChildren; ++i)
		grid.erase(m_pFaceChild[i]);
	clear();
}

void MGVolumeInfo::erase_all_children(Grid& grid)
{
	if(m_pVrtChild)
		grid.erase(m_pVrtChild);
	for(int i = 0; i < m_numEdgeChildren; ++i)
		grid.erase(m_pEdgeChild[i]);
	for(int i = 0; i < m_numFaceChildren; ++i)
		grid.erase(m_pFaceChild[i]);
	for(int i = 0; i < m_numVolChildren; ++i)
		grid.erase(m_pVolChild[i]);
	clear();
}

}//	end of namespace

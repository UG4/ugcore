/*
 * periodic_boundary_manager.cpp
 *
 *  Created on: 05.12.2012
 *      Author: marscher
 */
#include "periodic_boundary_manager.h"

namespace ug {

PeriodicBoundaryManager::~PeriodicBoundaryManager()
{
	for (VertexBaseIterator iter = m_pGrid->begin<VertexBase>();
			iter != m_pGrid->end<VertexBase>(); ++iter) {
		if (master(*iter)) delete m_aaGroupVRT[*iter];
	}

	for (EdgeBaseIterator iter = m_pGrid->begin<EdgeBase>();
			iter != m_pGrid->end<EdgeBase>(); ++iter) {
		if(master(*iter)) delete m_aaGroupEDG[*iter];
	}

	for (FaceIterator iter = m_pGrid->begin<Face>();
			iter != m_pGrid->end<Face>(); ++iter) {
		if(master(*iter)) delete m_aaGroupFCE[*iter];
	}

	for (VolumeIterator iter = m_pGrid->begin<Volume>();
			iter != m_pGrid->end<Volume>(); ++iter) {
		if(master(*iter)) delete m_aaGroupVOL[*iter];
	}

	// set grid to NULL to detach groups from grid
	set_grid(NULL);
}

void PeriodicBoundaryManager::set_grid(Grid* g)
{
	// group attachments
	Attachment<Group<VertexBase>*> aGroupVRT;
	Attachment<Group<EdgeBase>*> aGroupEDG;
	Attachment<Group<Face>*> aGroupFCE;
	Attachment<Group<Volume>*> aGroupVOL;

	if(g != NULL) {
		m_pGrid = g;

		m_pGrid->attach_to_vertices_dv(aGroupVRT, NULL);
		m_pGrid->attach_to_edges_dv(aGroupEDG, NULL);
		m_pGrid->attach_to_faces_dv(aGroupFCE, NULL);
		m_pGrid->attach_to_volumes_dv(aGroupVOL, NULL);

		// access grid with those attachments
		m_aaGroupVRT.access(*m_pGrid, aGroupVRT);
		m_aaGroupEDG.access(*m_pGrid, aGroupEDG);
		m_aaGroupFCE.access(*m_pGrid, aGroupFCE);
		m_aaGroupVOL.access(*m_pGrid, aGroupVOL);
	}

	// detach groups
	if(g == NULL) {
		m_pGrid->detach_from_vertices(aGroupVRT);
		m_pGrid->detach_from_edges(aGroupEDG);
		m_pGrid->detach_from_faces(aGroupFCE);
		m_pGrid->detach_from_volumes(aGroupVOL);
	}
}

template <>
const Grid::AttachmentAccessor<VertexBase, Attachment<PeriodicBoundaryManager::Group<VertexBase>* > >&
PeriodicBoundaryManager::get_group_accessor() const {
	return m_aaGroupVRT;
}

template <>
const Grid::AttachmentAccessor<EdgeBase, Attachment<PeriodicBoundaryManager::Group<EdgeBase>* > >&
PeriodicBoundaryManager::get_group_accessor() const {
	return m_aaGroupEDG;
}

template <>
const Grid::AttachmentAccessor<Face, Attachment<PeriodicBoundaryManager::Group<Face>* > >&
PeriodicBoundaryManager::get_group_accessor() const {
	return m_aaGroupFCE;
}

template <>
const Grid::AttachmentAccessor<Volume, Attachment<PeriodicBoundaryManager::Group<Volume>* > >&
PeriodicBoundaryManager::get_group_accessor() const {
	return m_aaGroupVOL;
}


} // end of namespace ug

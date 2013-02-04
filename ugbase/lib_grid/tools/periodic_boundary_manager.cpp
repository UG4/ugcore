/*
 * periodic_boundary_manager.cpp
 *
 *  Created on: 05.12.2012
 *      Author: marscher
 */
#include "periodic_boundary_manager.h"

namespace ug {

PeriodicBoundaryManager::PeriodicBoundaryManager() :
				m_pGrid(NULL), m_pSH(NULL) {}

PeriodicBoundaryManager::~PeriodicBoundaryManager()
{
	// set grid to NULL to detach groups from grid
	set_grid(NULL);
}

void PeriodicBoundaryManager::set_grid(Grid* g)
{
	// group attachments
	Attachment<Group<VertexBase>*> aGroupVRT;
	Attachment<Group<EdgeBase>*> aGroupEDG;
	Attachment<Group<Face>*> aGroupFCE;

	// periodic info attachments
	Attachment<PeriodicStatus> aPeriodicStatus;

	// detach groups and unregister observer
	if(m_pGrid) {
		m_pGrid->detach_from_vertices(aGroupVRT);
		m_pGrid->detach_from_edges(aGroupEDG);
		m_pGrid->detach_from_faces(aGroupFCE);

		m_pGrid->detach_from_vertices(aPeriodicStatus);
		m_pGrid->detach_from_edges(aPeriodicStatus);
		m_pGrid->detach_from_faces(aPeriodicStatus);

		m_pGrid->unregister_observer(this);
	}

	m_pGrid = dynamic_cast<MultiGrid*>(g);

	// attach groups and register observer
	if(m_pGrid != NULL) {
		m_pGrid->attach_to_vertices_dv(aGroupVRT, NULL);
		m_pGrid->attach_to_edges_dv(aGroupEDG, NULL);
		m_pGrid->attach_to_faces_dv(aGroupFCE, NULL);

		m_pGrid->attach_to_vertices_dv(aPeriodicStatus, P_NOT_PERIODIC);
		m_pGrid->attach_to_edges_dv(aPeriodicStatus, P_NOT_PERIODIC);
		m_pGrid->attach_to_faces_dv(aPeriodicStatus, P_NOT_PERIODIC);

		// access grid with those attachments
		m_aaGroupVRT.access(*m_pGrid, aGroupVRT);
		m_aaGroupEDG.access(*m_pGrid, aGroupEDG);
		m_aaGroupFCE.access(*m_pGrid, aGroupFCE);

		m_aaPeriodicStatusVRT.access(*m_pGrid, aPeriodicStatus);
		m_aaPeriodicStatusEDG.access(*m_pGrid, aPeriodicStatus);
		m_aaPeriodicStatusFCE.access(*m_pGrid, aPeriodicStatus);

		uint options = OT_GRID_OBSERVER | OT_VERTEX_OBSERVER | OT_EDGE_OBSERVER |
				OT_FACE_OBSERVER;
		m_pGrid->register_observer(this, options);
	}
}

Grid* PeriodicBoundaryManager::get_grid() const {
	return m_pGrid;
}

void PeriodicBoundaryManager::set_subset_handler(ISubsetHandler* sh) {
	m_pSH = sh;
	// very small memory overhead buys constant access time in match() methods
	m_vIdentifier.resize(m_pSH->num_subsets());
}

// group accessors
template <>
Grid::AttachmentAccessor<VertexBase, Attachment<PeriodicBoundaryManager::Group<VertexBase>* > >&
PeriodicBoundaryManager::get_group_accessor() {
	return m_aaGroupVRT;
}

template <>
Grid::AttachmentAccessor<EdgeBase, Attachment<PeriodicBoundaryManager::Group<EdgeBase>* > >&
PeriodicBoundaryManager::get_group_accessor() {
	return m_aaGroupEDG;
}

template <>
Grid::AttachmentAccessor<Face, Attachment<PeriodicBoundaryManager::Group<Face>* > >&
PeriodicBoundaryManager::get_group_accessor() {
	return m_aaGroupFCE;
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
	UG_THROW("not impled");
}

// periodic status accessors
template <>
const Grid::AttachmentAccessor<VertexBase, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	return m_aaPeriodicStatusVRT;
}

template <>
Grid::AttachmentAccessor<VertexBase, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() {
	return m_aaPeriodicStatusVRT;
}

template <>
const Grid::AttachmentAccessor<EdgeBase, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	return m_aaPeriodicStatusEDG;
}

template <>
Grid::AttachmentAccessor<EdgeBase, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() {
	return m_aaPeriodicStatusEDG;
}

template <>
Grid::AttachmentAccessor<Face, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() {
	return m_aaPeriodicStatusFCE;
}

template <>
const Grid::AttachmentAccessor<Face, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	return m_aaPeriodicStatusFCE;
}

template <>
Grid::AttachmentAccessor<Volume, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() {
	UG_THROW("not impled");
}

template <>
const Grid::AttachmentAccessor<Volume, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	UG_THROW("not impled");
}

///// grid observer methods
void PeriodicBoundaryManager::grid_to_be_destroyed(Grid* grid) {

	/// delete groups
	for (VertexBaseIterator iter = m_pGrid->begin<VertexBase>();
			iter != m_pGrid->end<VertexBase>(); ++iter) {
		if(is_master(*iter)) delete m_aaGroupVRT[*iter];
	}

	for (EdgeBaseIterator iter = m_pGrid->begin<EdgeBase>();
			iter != m_pGrid->end<EdgeBase>(); ++iter) {
		if(is_master(*iter)) delete m_aaGroupEDG[*iter];
	}

	for (FaceIterator iter = m_pGrid->begin<Face>();
			iter != m_pGrid->end<Face>(); ++iter) {
		if(is_master(*iter)) delete m_aaGroupFCE[*iter];
	}

	set_grid(NULL);
}

void PeriodicBoundaryManager::vertex_created(Grid* grid, VertexBase* vrt,
		GeometricObject* pParent, bool replacesParent) {
	handle_creation_cast_wrapper(vrt, pParent, replacesParent);
}

void PeriodicBoundaryManager::edge_created(Grid* grid, EdgeBase* e, GeometricObject* pParent,
		bool replacesParent) {
	handle_creation_cast_wrapper(e, pParent, replacesParent);
}

void PeriodicBoundaryManager::face_created(Grid* grid, Face* f, GeometricObject* pParent,
		bool replacesParent) {
	handle_creation_cast_wrapper(f, pParent, replacesParent);
}

void PeriodicBoundaryManager::vertex_to_be_erased(Grid* grid, VertexBase* vrt,
		VertexBase* replacedBy) {
	handle_deletion(vrt, replacedBy);
}

void PeriodicBoundaryManager::edge_to_be_erased(Grid* grid, EdgeBase* e,
		EdgeBase* replacedBy) {
	handle_deletion(e, replacedBy);
}

void PeriodicBoundaryManager::face_to_be_erased(Grid* grid, Face* f,
		Face* replacedBy) {
	handle_deletion(f, replacedBy);
}

void PeriodicBoundaryManager::set_identifier(IIdentifier*i, size_t si) {
	UG_ASSERT(m_vIdentifier.capacity() >= si, "identifier vector not big enough")
	m_vIdentifier[si] = SmartPtr<IIdentifier>(i);
}

} // end of namespace ug

/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Martin Scherer
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

#include "periodic_boundary_manager.h"

namespace ug {


PeriodicBoundaryManager::PeriodicBoundaryManager() :
				m_pGrid(nullptr)/*, m_pSH(nullptr)*/ {}

PeriodicBoundaryManager::~PeriodicBoundaryManager()
{
	// set grid to nullptr to detach groups from grid
	set_grid(nullptr);
}

void PeriodicBoundaryManager::set_grid(Grid* g)
{
	// group attachments
	Attachment<Group<Vertex>*> aGroupVRT;
	Attachment<Group<Edge>*> aGroupEDG;
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
	if(m_pGrid != nullptr) {
		m_pGrid->attach_to_vertices_dv(aGroupVRT, nullptr);
		m_pGrid->attach_to_edges_dv(aGroupEDG, nullptr);
		m_pGrid->attach_to_faces_dv(aGroupFCE, nullptr);

		m_pGrid->attach_to_vertices_dv(aPeriodicStatus, PeriodicStatus::P_NOT_PERIODIC);
		m_pGrid->attach_to_edges_dv(aPeriodicStatus, PeriodicStatus::P_NOT_PERIODIC);
		m_pGrid->attach_to_faces_dv(aPeriodicStatus, PeriodicStatus::P_NOT_PERIODIC);

		// access grid with those attachments
		m_aaGroupVRT.access(*m_pGrid, aGroupVRT);
		m_aaGroupEDG.access(*m_pGrid, aGroupEDG);
		m_aaGroupFCE.access(*m_pGrid, aGroupFCE);

		m_aaPeriodicStatusVRT.access(*m_pGrid, aPeriodicStatus);
		m_aaPeriodicStatusEDG.access(*m_pGrid, aPeriodicStatus);
		m_aaPeriodicStatusFCE.access(*m_pGrid, aPeriodicStatus);

		ObserverType options = ObserverType::OT_GRID_OBSERVER | ObserverType::OT_VERTEX_OBSERVER | ObserverType::OT_EDGE_OBSERVER |
				ObserverType::OT_FACE_OBSERVER;
		m_pGrid->register_observer(this, options);
	}
}

Grid* PeriodicBoundaryManager::get_grid() const {
	return m_pGrid;
}

//void PeriodicBoundaryManager::set_subset_handler(ISubsetHandler* sh) {
//	m_pSH = sh;
//	// very small memory overhead buys constant access time in match() methods
//	m_vIdentifier.resize(m_pSH->num_subsets());
//}

// group accessors
template <>
Grid::AttachmentAccessor<Vertex, Attachment<PeriodicBoundaryManager::Group<Vertex>* > >&
PeriodicBoundaryManager::get_group_accessor() {
	return m_aaGroupVRT;
}

template <>
Grid::AttachmentAccessor<Edge, Attachment<PeriodicBoundaryManager::Group<Edge>* > >&
PeriodicBoundaryManager::get_group_accessor() {
	return m_aaGroupEDG;
}

template <>
Grid::AttachmentAccessor<Face, Attachment<PeriodicBoundaryManager::Group<Face>* > >&
PeriodicBoundaryManager::get_group_accessor() {
	return m_aaGroupFCE;
}

template <>
const Grid::AttachmentAccessor<Vertex, Attachment<PeriodicBoundaryManager::Group<Vertex>* > >&
PeriodicBoundaryManager::get_group_accessor() const {
	return m_aaGroupVRT;
}

template <>
const Grid::AttachmentAccessor<Edge, Attachment<PeriodicBoundaryManager::Group<Edge>* > >&
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
	UG_THROW("not implemented");
}

// periodic status accessors
template <>
const Grid::AttachmentAccessor<Vertex, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	return m_aaPeriodicStatusVRT;
}

template <>
Grid::AttachmentAccessor<Vertex, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() {
	return m_aaPeriodicStatusVRT;
}

template <>
const Grid::AttachmentAccessor<Edge, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	return m_aaPeriodicStatusEDG;
}

template <>
Grid::AttachmentAccessor<Edge, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
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
	UG_THROW("not implemented");
}

template <>
const Grid::AttachmentAccessor<Volume, Attachment<PeriodicBoundaryManager::PeriodicStatus> >&
PeriodicBoundaryManager::get_periodic_status_accessor() const {
	UG_THROW("not implemented");
}

///// grid observer methods
void PeriodicBoundaryManager::grid_to_be_destroyed(Grid* grid) {

	/// delete groups
	for (VertexIterator iter = m_pGrid->begin<Vertex>();
			iter != m_pGrid->end<Vertex>(); ++iter) {
		if(is_master(*iter)) delete m_aaGroupVRT[*iter];
	}

	for (EdgeIterator iter = m_pGrid->begin<Edge>();
			iter != m_pGrid->end<Edge>(); ++iter) {
		if(is_master(*iter)) delete m_aaGroupEDG[*iter];
	}

	for (FaceIterator iter = m_pGrid->begin<Face>();
			iter != m_pGrid->end<Face>(); ++iter) {
		if(is_master(*iter)) delete m_aaGroupFCE[*iter];
	}

	set_grid(nullptr);
}

void PeriodicBoundaryManager::vertex_created(Grid* grid, Vertex* vrt,
		GridObject* pParent, bool replacesParent) {
	handle_creation_cast_wrapper(vrt, pParent, replacesParent);
}

void PeriodicBoundaryManager::edge_created(Grid* grid, Edge* e, GridObject* pParent,
		bool replacesParent) {
	handle_creation_cast_wrapper(e, pParent, replacesParent);
}

void PeriodicBoundaryManager::face_created(Grid* grid, Face* f, GridObject* pParent,
		bool replacesParent) {
	handle_creation_cast_wrapper(f, pParent, replacesParent);
}

void PeriodicBoundaryManager::vertex_to_be_erased(Grid* grid, Vertex* vrt,
		Vertex* replacedBy) {
	handle_deletion(vrt, replacedBy);
}

void PeriodicBoundaryManager::edge_to_be_erased(Grid* grid, Edge* e,
		Edge* replacedBy) {
	handle_deletion(e, replacedBy);
}

void PeriodicBoundaryManager::face_to_be_erased(Grid* grid, Face* f,
		Face* replacedBy) {
	handle_deletion(f, replacedBy);
}

/**
 * performs following checks on all elements in goc1 and goc2
 * 1. all elements given are periodic
 * 2. masters of groups are valid (master pointer of group valid + have at least one child)
 * 3. no duplicate master, slave pairs exists
 * 4.
 */
bool PeriodicBoundaryManager::check_periodicity(
		const GridObjectCollection& goc1,
		const GridObjectCollection& goc2, ISubsetHandler* sh) {

	Group<Vertex>::unique_pairs s_vert;
	Group<Edge>::unique_pairs s_edge;
	Group<Face>::unique_pairs s_face;

	UG_ASSERT(goc1.num_levels() == goc2.num_levels(),
			"collections have different mg levels!")

	for (size_t lvl = 0; lvl < goc1.num_levels(); lvl++) {
		// check faces
		check_elements_periodicity<Face>(goc1.begin<Face>(lvl),
				goc1.end<Face>(lvl), s_face, sh);
		check_elements_periodicity<Face>(goc2.begin<Face>(lvl),
				goc2.end<Face>(lvl), s_face, sh);

		// check edges
		check_elements_periodicity<Edge>(goc1.begin<Edge>(lvl),
				goc1.end<Edge>(lvl), s_edge, sh);
		check_elements_periodicity<Edge>(goc2.begin<Edge>(lvl),
				goc2.end<Edge>(lvl), s_edge, sh);

		// check vertices
		check_elements_periodicity<Vertex>(goc1.begin<Vertex>(lvl),
				goc1.end<Vertex>(lvl), s_vert, sh);
		check_elements_periodicity<Vertex>(goc2.begin<Vertex>(lvl),
				goc2.end<Vertex>(lvl), s_vert, sh);
	}

	return true;
}


template <typename elem_t>
void PeriodicBoundaryManager::validity_check()
{
	using iter_t = typename Grid::traits<elem_t>::iterator;
	using slave_container_t = typename Group<elem_t>::SlaveContainer;
	using slave_iter_t = typename slave_container_t::iterator;

	Grid::AttachmentAccessor<elem_t, Attachment<PeriodicStatus> >
		aaPS = get_periodic_status_accessor<elem_t>();

	Grid::AttachmentAccessor<elem_t, Attachment<Group<elem_t>* > >
		aaGroup = get_group_accessor<elem_t>();

	Grid& g = *m_pGrid;

	for(iter_t i_elem = g.begin<elem_t>(); i_elem != g.end<elem_t>(); ++ i_elem)
	{
		elem_t* e = *i_elem;
		if(aaGroup[e]){
			bool isMaster = (aaGroup[e]->m_master == e);
			int elemSlaveCount = 0;
			int totalSlaveCount = 0;
			slave_container_t* slaveCon = &aaGroup[e]->get_slaves();
			for(slave_iter_t i = slaveCon->begin(); i != slaveCon->end(); ++i)
			{
				++totalSlaveCount;
				if(*i == e)
					++elemSlaveCount;
			}

			if(aaGroup[e]->m_master == nullptr){
				UG_THROW("Group doesn't contain a master!\n"
						 << "  elem-info: " << ElementDebugInfo(g, e));
			}
			if(totalSlaveCount == 0){
				UG_THROW("Group doesn't contain slaves!\n"
						 << "  elem-info: " << ElementDebugInfo(g, e));
			}

			if(isMaster && (elemSlaveCount > 0)){
				UG_THROW("An element is contained both in the master list and in the"
						 " slave list of its group!\n"
						 << "  elem-info: " << ElementDebugInfo(g, e));
			}
			if(elemSlaveCount > 1){
				UG_THROW("Multi-occurrence of an element in a groups slave list!\n"
						 << "  elem-info: " << ElementDebugInfo(g, e));
			}
			if((!isMaster) && (elemSlaveCount == 0)){
				UG_THROW("An element points to a group but is not contained in that group!\n"
						 << "  elem-info: " << ElementDebugInfo(g, e));
			}
			
		//	the element is either master or slave. now check whether the states are correct
			if(isMaster){
				if(aaPS[e] != P_MASTER){
					UG_THROW("An element is a groups master, but its state does not say so!\n"
						 	 << "  state: " << aaPS[e] << "\n"
					 		 << "  elem-info: " << ElementDebugInfo(g, e));
				}
			//	make sure that all slaves have their group assigned
				for(slave_iter_t i = slaveCon->begin(); i != slaveCon->end(); ++i)
				{
					if(aaGroup[*i] != aaGroup[e]){
						UG_THROW("A slave element was encountered, which does not"
								 << " point to the same group as it's master!\n"
								 << "  master: " << ElementDebugInfo(g, e) << "\n"
								 << "  slave: " << ElementDebugInfo(g, *i) << "\n");
					}
				}
			}
			else{
				if(aaPS[e] == P_SLAVE_MASTER_UNKNOWN){
					UG_THROW("P_SLAVE_MASTER_UNKNOWN is not a valid state!\n"
					 		 << "  elem-info: " << ElementDebugInfo(g, e));
				}
				else if(aaPS[e] != P_SLAVE)
				{
					UG_THROW("An element is a slave, but its state does not say so!\n"
						 	 << "  state: " << aaPS[e] << "\n"
					 		 << "  elem-info: " << ElementDebugInfo(g, e));
				}

				if(aaGroup[e] != aaGroup[aaGroup[e]->m_master]){
					UG_THROW("A master element was encountered, which does not"
							 << " point to the same group as it's slave!\n"
							 << "  master: " << ElementDebugInfo(g, aaGroup[e]->m_master) << "\n"
							 << "  slave: " << ElementDebugInfo(g, e) << "\n");
				}
			}
		}
		else if(aaPS[e] != P_NOT_PERIODIC){
			UG_THROW("An element is marked as periodic but not contained in a group!\n"
					 << "  state: " << aaPS[e] << "\n"
					 << "  elem-info: " << ElementDebugInfo(g, e));
		}
	}
}

void PeriodicBoundaryManager::validity_check()
{
	validity_check<Vertex>();
	validity_check<Edge>();
	validity_check<Face>();
}

//void PeriodicBoundaryManager::set_identifier(SmartPtr<IIdentifier> i, size_t si) {
//	UG_ASSERT(m_vIdentifier.capacity() >= si, "identifier vector not big enough")
//	m_vIdentifier[si] = i;
//}

} // end of namespace ug

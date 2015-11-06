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

#ifndef PERIODIC_IDENTIFIER_IMPL_HPP_
#define PERIODIC_IDENTIFIER_IMPL_HPP_

// include declarations
#include "periodic_boundary_manager.h"
#include "lib_disc/domain.h"
#include "lib_grid/algorithms/debug_util.h"
#include "common/assert.h"
#include "common/error.h"
#include "pcl/pcl_base.h"

#include <boost/mpl/map.hpp>
#include <boost/mpl/at.hpp>

#include <algorithm>

namespace ug {

template <class TAAPos>
template <class TElem>
bool ParallelShiftIdentifier<TAAPos>::match_impl(TElem* e1, TElem* e2) const {
	if (e1 == e2)
		return false;

	AttachmentType c1 = CalculateCenter(e1, m_aaPos),
			c2 = CalculateCenter(e2, m_aaPos), diff, error;
	bool result = false;

	VecSubtract(diff, c1, c2);
	VecSubtract(error, diff, m_shift);
	number len = VecLengthSq(error);
	if (std::abs(len) < 10E-8)
		result = true;
	else // check for opposite shift
	{
		VecSubtract(error, diff, m_shift_opposite);
		len = VecLengthSq(error);
		if (std::abs(len) < 10E-8)
			result = true;
	}

	return result;
}

template <class TElem>
void PeriodicBoundaryManager::identify(TElem* e1, TElem* e2,
		IIdentifier& ident) {
	typedef typename
			Grid::traits<typename TElem::side>::secure_container container;

	// determine masters
	TElem* m1 = master(e1);
	TElem* m2 = master(e2);

	if (m1 || m2) {
		//	if m1 == m2, there's nothing to do
		if (m1 != m2) {
			if (!m1) { // m2 is master
				make_slave(group(m2), e1);
			} else if (!m2) { // m1 is master
				make_slave(group(m1), e2);
			} else if (m1 && m2) { // both are distinct masters
				merge_groups(group(m1), group(m2));
			} else {
				UG_THROW("should never get here")
			}
		}
	} else {
		// create new group with e1 as master
		Group<TElem>* g = new Group<TElem>(e1);
		// set group of master
		set_group(g, e1);
		// make e2 slave
		make_slave(g, e2);
	}

	// while elements have sides, recursively identify sub type elements
	if (TElem::HAS_SIDES) {
		container sides1, sides2;
		// collect sides and identify them
		m_pGrid->associated_elements<TElem>(sides1, e1);
		m_pGrid->associated_elements<TElem>(sides2, e2);
		for (size_t i = 0; i < sides1.size(); ++i) {
			for (size_t j = 0; j < sides2.size(); ++j) {
				if(ident.match(sides1[i], sides2[j])) {
					identify(sides1[i], sides2[j], ident);
				}
			}
		}
	}
}

// is considered periodic if it is a master or a slave
template <class TElem>
bool PeriodicBoundaryManager::is_periodic(TElem* e) const {
	if(e)
		return get_periodic_status_accessor<TElem>()[e] != P_NOT_PERIODIC;
	else
		UG_THROW("null pointer is never periodic.");
}

// gets master of e, may be null
template <class TElem>
TElem* PeriodicBoundaryManager::master(TElem* e) const {
	if (group(e)) {
		return group(e)->m_master;
	}
	return NULL;
}

// gets slaves of e
template <class TElem>
typename PeriodicBoundaryManager::Group<TElem>::SlaveContainer*
PeriodicBoundaryManager::slaves(TElem* e) const {
	if (group(e)) {
		return &group(e)->get_slaves();
	}
	return NULL;
}

template <class TElem>
void PeriodicBoundaryManager::print_identification() const {
	typedef typename ElementStorage<TElem>::SectionContainer::iterator Iterator;
	typedef typename Group<TElem>::SlaveIterator SlaveIter;

	for (Iterator elem = m_pGrid->begin<TElem>(); elem != m_pGrid->end<TElem>();
			++elem) {
		// log masters and their slaves
		if (!(master(*elem) == *elem))
			continue;

		Group<TElem>* g = group(*elem);
		UG_ASSERT(g, "group not valid")
		UG_LOG("group of " << (*elem)->reference_object_id() << "\tlevel: " <<
				m_pGrid->get_level(*elem) << "\tmaster: " <<
				GetGridObjectCenter(*m_pGrid, g->m_master) << "\tslaves: ");
		for (SlaveIter slave = g->get_slaves().begin();
				slave != g->get_slaves().end(); ++slave) {
			TElem* e = *slave;
			UG_LOG(GetGridObjectCenter(*m_pGrid, e) << ", ")
			UG_ASSERT(m_pGrid->get_level(*elem) == m_pGrid->get_level(e),
					"wrong level in group")
		}
		UG_LOG(std::endl)
	}
}

/**
 * TParent should be only of type Vertex, Edge, Face.
 * Volumes are not meant to be periodic.
 *
 * If replacesParent is true, e is meant to replace pParent
 */
template <class TElem>
void PeriodicBoundaryManager::replace_parent(TElem* e, TElem* pParent) {

	if (is_master(pParent)) {
		// if parent is master, set newly created item as new master
		Group<TElem>* g = group(pParent);
		g->m_master = e;
		set_group(g, e);
	} else if(is_slave(pParent)) { // slave
	//	iterate over parent-groups slave list and replace pointer to parent
		typename Group<TElem>::SlaveContainer* slaveCon = slaves(pParent);
		for(typename Group<TElem>::SlaveContainer::iterator i = slaveCon->begin();
			i != slaveCon->end(); ++i)
		{
			if(*i == pParent)
				*i = e;
		}
		set_group(group(pParent), e);
	}
}

template <class TElem, class TParent>
void PeriodicBoundaryManager::handle_creation(TElem* e, TParent* pParent) {

	typedef typename Group<TParent>::SlaveContainer ParentSlaveContainer;

	if (!pParent)
		return;

	if (!is_periodic(pParent))
		return;

	UG_ASSERT(m_pGrid,
			"PeriodicBoundaryManager has to operate on a grid. No grid assigned.");
	MultiGrid& mg = *m_pGrid;

	mg.begin_marking();
	if(TElem::BASE_OBJECT_ID != VERTEX){
		Grid::vertex_traits::secure_container vrts;
		mg.associated_elements(vrts, e);
		for(size_t i = 0; i < vrts.size(); ++i){
			Group<Vertex>* grp = group(vrts[i]);
			if(!grp)
				continue;
			mg.mark(grp->m_master);
			mg.mark(grp->get_slaves().begin(), grp->get_slaves().end());
		}
	}

	if (is_master<TParent>(pParent)) {
		// create new group for e, with e as master
		Group<TElem>* newGroup = new Group<TElem>(e);
		UG_ASSERT(group(e) == NULL, "element already has a group.")
		set_group(newGroup, e);

		// iterate over slaves of parent
		ParentSlaveContainer* parentSlaves = slaves(pParent);
		for (typename ParentSlaveContainer::iterator i_slave =
				parentSlaves->begin(); i_slave != parentSlaves->end();
				++i_slave) {
			TParent* parentSlave = *i_slave;
			UG_ASSERT(parentSlave, "parent slave not valid")

			// iterate over all children of parentSlave and make them slaves of e
			// only matching children of type TElem are considered.
			// note that not all children already have to exist at this point.
			// children which are created later on are handled by the 'else' section.
			for (size_t i_child = 0;
					i_child < mg.num_children<TElem>(parentSlave); ++i_child) {
				TElem* c = mg.get_child<TElem>(parentSlave, i_child);
				UG_ASSERT(!(e->base_object_id() == VERTEX)
						  || (mg.num_children<TElem>(parentSlave) == 1),
						  "At most one 1 vertex-child is currently allowed");
				// We use a special case for vertices here, since the position
				// attachment has not yet been set (this is currently done after
				// all observers have been informed about the creation of the new
				// vertex). Note that this is no problem for edges or faces, since
				// associated vertices have been properly initialized at that time.
				if((e->base_object_id() == VERTEX)) {
					make_slave(newGroup, c);
				} else {
					Grid::vertex_traits::secure_container vrts;
					mg.associated_elements(vrts, c);
					bool allMarked = true;
					for(size_t i = 0; i < vrts.size(); ++i){
						if(!mg.is_marked(vrts[i])){
							allMarked = false;
							break;
						}
					}

					if(allMarked)
						make_slave(newGroup, c);
				}
			}
		}
	} else {
		// find the associated master of e by checking for a matching child
		// in the set of children of the parents master.
		// If no appropriate child already exists at this point, no master can
		// be found. We simply leave the new element as it is, since it will
		// be found as a slave when the associated master will be created later
		// on (the 'if(is_master...)' case applies then.

		// attention: parentGroup may be null
		Group<TParent>* parentGroup = group<TParent>(pParent);

		if (parentGroup == NULL) {
			get_periodic_status_accessor<TElem>()[e] = P_SLAVE_MASTER_UNKNOWN;
			return;
		}

		TParent* parentMaster = parentGroup->m_master;
		bool master_found = false;
		for (size_t i_child = 0; i_child < mg.num_children<TElem>(parentMaster);
				++i_child) {
			TElem* c = mg.get_child<TElem>(parentMaster, i_child);
			UG_ASSERT(!(e->base_object_id() == VERTEX)
					  || (mg.num_children<TElem>(parentMaster) == 1),
					  "At most one 1 vertex-child is currently allowed");
			// We use a special case for vertices here, since the position
			// attachment has not yet been set (this is currently done after
			// all observers have been informed about the creation of the new
			// vertex). Note that this is no problem for edges or faces, since
			// associated vertices have been properly initialized at that time.
			if((e->base_object_id() == VERTEX)) {
				make_slave(group(c), e);
				master_found = true;
				break;
			} else {
				Grid::vertex_traits::secure_container vrts;
				mg.associated_elements(vrts, c);
				bool allMarked = true;
				for(size_t i = 0; i < vrts.size(); ++i){
					if(!mg.is_marked(vrts[i])){
						allMarked = false;
						break;
					}
				}

				if(allMarked){
					make_slave(group(c), e);
					master_found = true;
					break;
				}
			}
		}

		// wait until a master for e will be known
		if (!master_found) {
			get_periodic_status_accessor<TElem>()[e] = P_SLAVE_MASTER_UNKNOWN;
		}
	}

	mg.end_marking();
}

template <class TElem>
void PeriodicBoundaryManager::handle_creation_cast_wrapper(TElem* e,
		GridObject* pParent, bool replacesParent) {
	// we can only identify periodic elements, which have a periodic parent
	if(!pParent)
		return;

	if(replacesParent){
		replace_parent(e, static_cast<TElem*>(pParent));
	}
	else{
		switch (pParent->base_object_id()) {
		case VERTEX:
			handle_creation(e, static_cast<Vertex*>(pParent));
			break;
		case EDGE:
			handle_creation(e, static_cast<Edge*>(pParent));
			break;
		case FACE:
			handle_creation(e, static_cast<Face*>(pParent));
			break;
		// ignore volumes, as these are not meant to be periodic
		case VOLUME:
			break;
		default:
			UG_THROW("no handling for parent type: " << pParent->base_object_id())
		}
	}
}

/// handles deletion of element type
template <class TElem>
void PeriodicBoundaryManager::handle_deletion(TElem* e, TElem* replacedBy) {
	if((!is_periodic(e)) || replacedBy)
		return;

	if (is_master(e)) {
		// delete a master completely...
		if (replacedBy) {
			UG_ASSERT(!group(replacedBy),
			"replacing element is already in group")
			make_master(group(e), replacedBy);
		} else {
			remove_group(group(e));
		}
	} else { // slave
		UG_ASSERT(is_slave(e), "e should be a slave")
		if(!remove_slave(e))
			UG_THROW("old slave not removed.")
	}
}

template <class TElem>
void PeriodicBoundaryManager::make_slave(Group<TElem>* g, TElem* slave) {
	UG_ASSERT(g, "invalid group")
	UG_ASSERT(slave, "invalid slave")
	UG_ASSERT(group(slave) != g, "trying to add a duplicate slave!")
	UG_ASSERT(m_pGrid->get_level(g->m_master) == m_pGrid->get_level(slave),
			"level of slave and group mismatch")

	// if slave is already engrouped, remove it first
	if(group(slave)) {
		UG_LOG("slave already engrouped. removing it from former group\n")
		if(!remove_slave(slave)) {
			UG_THROW("slave could not be removed")
		}
	}

	// add slave to group and set group attachment/status
	g->add_slave(slave);
	set_group(g, slave);
}

template <class TElem>
void PeriodicBoundaryManager::make_master(Group<TElem>* g, TElem* new_master) {
	// group has a master, reset its group
	if(g->m_master) {
		set_group<TElem>(NULL, g->m_master);
	}
	g->m_master = new_master;
}

template <class TElem>
bool PeriodicBoundaryManager::remove_slave(TElem* e) {
	if (slaves(e)) {
		typename Group<TElem>::SlaveContainer& s = *slaves(e);
		typename Group<TElem>::SlaveIterator pos = std::find(s.begin(), s.end(), e);
		if (pos != s.end()) {
//			s.erase(pos);
//			set_group<TElem>(NULL, e);

			// todo this is the only slave left, remove whole group?!
			if(s.size() == 1) {
				UG_LOGN("delete last slave, remove_group")
				remove_group(group(e));
			}
			// todo what if e was the last slave of of group?
//			if(s.empty()) {
//				UG_THROW("empty group leftover....")
//			}
			return true;
		}
	}
	return false;
}

template <class TElem>
void PeriodicBoundaryManager::remove_group(Group<TElem>* g) {
	UG_ASSERT(g, "should remove invalid group")

	// reset group pointers of all group members to NULL
	set_group<TElem>(NULL, g->m_master);
	typename Group<TElem>::SlaveContainer& s = g->get_slaves();
	typename Group<TElem>::SlaveIterator iter;
	for (iter = s.begin(); iter != s.end(); ++iter) {
		set_group<TElem>(NULL, *iter);
	}

	delete g;
}

/**
 * merges g1 in g0 and deletes g1 afterwards
 */
template <class TElem>
void PeriodicBoundaryManager::merge_groups(Group<TElem>* g0, Group<TElem>* g1) {
	UG_ASSERT(g0 && g1, "groups not valid")
	UG_ASSERT(g0 != g1, "groups are equal")

	typedef typename Group<TElem>::SlaveContainer SlaveContainer;
	typedef typename Group<TElem>::SlaveIterator SlaveIterator;

	SlaveContainer& slaves_g1 = g1->get_slaves();

	// insert slaves of g1 at the end of slaves of g0
	for (SlaveIterator iter = slaves_g1.begin(); iter != slaves_g1.end();
			++iter) {
		TElem* e = *iter;
		UG_ASSERT(e, "slave not valid")
		UG_ASSERT(e != g0->m_master, "slave of g1 is master of g0!")
		// we do not use make_slave() here, because g1 will be deleted at the end
		g0->add_slave(e);
		set_group(g0, e);
	}

	// make old master a slave of group g1
	// we do not use make_slave() here, because g1 will be deleted at the end
	g0->add_slave(g1->m_master);
	set_group(g0, g1->m_master);

	// remove old group
	delete g1;
}

template <class TElem>
bool PeriodicBoundaryManager::is_slave(TElem* e) const {
	PeriodicStatus p = get_periodic_status_accessor<TElem>()[e];
	return (p == P_SLAVE || p == P_SLAVE_MASTER_UNKNOWN);
}

template <class TElem>
bool PeriodicBoundaryManager::is_master(TElem* e) const {
	return get_periodic_status_accessor<TElem>()[e] == P_MASTER;
}

template <class TElem>
PeriodicBoundaryManager::Group<TElem>* PeriodicBoundaryManager::group(
		TElem* e) const {
	UG_ASSERT(e, "element not valid.")
	return get_group_accessor<TElem>()[e];
}

template <class TElem>
void PeriodicBoundaryManager::set_group(Group<TElem>* g, TElem* e) {
	UG_ASSERT(e, "element not valid for attachment access.")
	// set group pointer of element e
	get_group_accessor<TElem>()[e] = g;

	// set periodic status
	if (g == NULL) {
		get_periodic_status_accessor<TElem>()[e] = P_NOT_PERIODIC;
	} else {
		if (g->m_master == e)
			get_periodic_status_accessor<TElem>()[e] = P_MASTER;
		else
			get_periodic_status_accessor<TElem>()[e] = P_SLAVE;
	}
}

template <class TElem, class TIterator>
void PeriodicBoundaryManager::check_elements_periodicity(
		TIterator begin,
		TIterator end,
		typename Group<TElem>::unique_pairs& s,
		ISubsetHandler* sh) {

	typedef typename Group<TElem>::SlaveContainer Container;
	typedef typename Group<TElem>::SlaveIterator SlaveIter;
	typedef typename ElementStorage<TElem>::SectionContainer::const_iterator SecContainerIter;

	for (SecContainerIter iter = begin; iter != end; ++iter) {
		TElem* e = *iter;
		if(! e)
			continue;
		if(!is_periodic(e)) {
			// lookup subset name of element
			const char* sh_name = "";
			if(sh) {
				int element_si = sh->get_subset_index(e);
				sh_name = sh->get_subset_name(element_si);
			}
			UG_THROW("Element in subset '" << sh_name
					<< "' is not periodic after identification: "
					<< GetGridObjectCenter(*get_grid(), e)
					<< "\nCheck your geometry for symmetry!\n")
		}

		if (master(e) == e) {
			Container* _slaves = slaves(e);
			if(! _slaves)
				UG_THROW("masters slave storage is not valid.")

			if(_slaves->empty())
				UG_THROW("master has no slaves")

			for (SlaveIter i = _slaves->begin(); i != _slaves->end(); ++i) {
				TElem* slave = *i;
				typename Group<TElem>::master_slave_pair p = std::make_pair(e, slave);
				bool inserted = (s.insert(p)).second;
				if(! inserted)
					UG_THROW("master/slave pair already exists.");
			}
		}
	}
}

template <class TDomain>
void IdentifySubsets(TDomain& dom, const char* sName1, const char* sName2) {
	// get subset handler from domain
	typedef typename TDomain::subset_handler_type subset_handler_type;

	subset_handler_type& sh = *dom.subset_handler();

	int si1 = sh.get_subset_index(sName1);
	int si2 = sh.get_subset_index(sName2);

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
		UG_THROW("sorry, in real parallel environment periodic bnds are not impled yet.")
#endif

	if (si1 == -1)
		UG_THROW("IdentifySubsets: given subset name " << sName1 << " does not exist");
	if (si2 == -1)
		UG_THROW("IdentifySubsets: given subset name " << sName2 << " does not exist");

	IdentifySubsets(dom, si1, si2);
}

/// performs geometric ident of periodic elements and master slave
template <class TDomain>
void IdentifySubsets(TDomain& dom, int sInd1, int sInd2) {

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1)
		UG_THROW("sorry, in real parallel environment periodic bnds are not impled yet.")
#endif

	if (sInd1 == -1 || sInd2 == -1) {
		UG_THROW("IdentifySubsets: at least one invalid subset given!")
	}

	if(sInd1 == sInd2) {
		UG_THROW("IdentifySubsets: can not identify two identical subsets!")
	}

	// ensure grid has support for periodic boundaries
	if (!dom.grid()->has_periodic_boundaries()) {
		dom.grid()->set_periodic_boundaries(true);
	}

	PeriodicBoundaryManager& pbm = *dom.grid()->periodic_boundary_manager();

	typedef typename TDomain::position_type position_type;
	typedef typename TDomain::position_accessor_type position_accessor_type;

	// get subset handler from domain
	typedef typename TDomain::subset_handler_type subset_handler_type;

	subset_handler_type& sh = *dom.subset_handler();

	// get aaPos from domain
	position_accessor_type& aaPos = dom.position_accessor();

	// create parallel shift identifier to match subset elements
	ParallelShiftIdentifier<position_accessor_type> ident(aaPos);

	// shift vector between subsets
	position_type shift;

	// collect all geometric objects (even for all levels in case of multi grid)
	GridObjectCollection goc1 = sh.get_grid_objects_in_subset(sInd1);
	GridObjectCollection goc2 = sh.get_grid_objects_in_subset(sInd2);

	if(goc1.num<Vertex>() != goc2.num<Vertex>()) {
		UG_THROW("IdentifySubsets: Given subsets have different number of vertices."
				"\nnum# in " << sh.get_subset_name(sInd1) << ": " << goc1.num<Vertex>() <<
				"\nnum# in " << sh.get_subset_name(sInd2) << ": " << goc2.num<Vertex>())
	}

	if(goc1.num<Edge>() != goc2.num<Edge>()) {
		UG_THROW("IdentifySubsets: Given subsets have different number of edges."
						"\nnum# in " << sh.get_subset_name(sInd1) << ": " << goc1.num<Edge>() <<
						"\nnum# in " << sh.get_subset_name(sInd2) << ": " << goc2.num<Edge>())
	}

	if(goc1.num<Face>() != goc2.num<Face>()) {
		UG_THROW("IdentifySubsets: Given subsets have different number of faces."
						"\nnum# in " << sh.get_subset_name(sInd1) << ": " << goc1.num<Face>() <<
						"\nnum# in " << sh.get_subset_name(sInd2) << ": " << goc2.num<Face>())
	}

	// map start type of recursion dependent to TDomain
	// in 3d start with faces, in 2d with edges, in 1d with vertices
	namespace mpl = boost::mpl;
	typedef		mpl::map<mpl::pair<Domain1d, Vertex>,
						 mpl::pair<Domain2d, Edge>,
						 mpl::pair<Domain3d, Face> > m;

	typedef typename mpl::at<m, TDomain>::type TElem;
	typedef typename ElementStorage<TElem>::SectionContainer::iterator gocIter;

	// calculate shift vector for top level
	position_type c1 = CalculateCenter(goc1.begin<TElem>(0), goc1.end<TElem>(0),
			aaPos);
	position_type c2 = CalculateCenter(goc2.begin<TElem>(0), goc2.end<TElem>(0),
			aaPos);

	VecSubtract(shift, c1, c2);
	ident.set_shift(shift);

	// for each level of multi grid. In case of simple grid only one iteration
	for (size_t lvl = 0; lvl < goc1.num_levels(); lvl++) {
		// identify corresponding elements for second subset. A element is considered
		// to have symmetric element in second subset if there exists a shift vector between them.
		// todo use kd-tree for fast lookup of shifted elements
		for (gocIter iter1 = goc1.begin<TElem>(lvl);
				iter1 != goc1.end<TElem>(lvl); ++iter1)
			for (gocIter iter2 = goc2.begin<TElem>(lvl);
					iter2 != goc2.end<TElem>(lvl); ++iter2) {
				if(ident.match(*iter1, *iter2)) {
					pbm.identify(*iter1, *iter2, ident);
				}
			}
	}

	// ensure periodic identification has been performed correctly
	pbm.check_periodicity(goc1, goc2, &sh);
}

} // end namespace ug

#endif /* PERIODIC_IDENTIFIER_IMPL_HPP_ */

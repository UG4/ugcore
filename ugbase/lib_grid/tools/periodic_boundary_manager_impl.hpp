/*
 * periodic_identifier_impl.hpp
 *
 *  Created on: 26.11.2012
 *      Author: marscher
 */

#ifndef PERIODIC_IDENTIFIER_IMPL_HPP_
#define PERIODIC_IDENTIFIER_IMPL_HPP_

// include declarations
#include "./periodic_boundary_manager.h"
#include "lib_grid/subset_handler.h"
#include "lib_disc/domain.h"
#include "lib_grid/algorithms/debug_util.h"
#include "common/assert.h"

#include <boost/mpl/map.hpp>
#include <boost/mpl/at.hpp>

#include <algorithm>
#include <set>

namespace ug {

template <class TAAPos>
template <class TElem>
bool ParallelShiftIdentifier<TAAPos>::match_impl(TElem* e1, TElem* e2) const {
	if (e1 == e2)
		return false;

	AttachmentType c1 = CalculateCenter(e1, m_aaPos), c2 = CalculateCenter(e2,
			m_aaPos), diff, error;
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
void PeriodicBoundaryManager::identifiy(TElem* e1, TElem* e2,
		IIdentifier* identifier) {
	typedef typename Grid::traits<typename TElem::side>::secure_container container;

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
	if (identifier != NULL || m_pIdentifier.get() != NULL) {
		if (TElem::HAS_SIDES) {
			container sides1, sides2;
			// collect sides and identify them
			m_pGrid->associated_elements<TElem>(sides1, e1);
			m_pGrid->associated_elements<TElem>(sides2, e2);
			for (size_t i = 0; i < sides1.size(); ++i) {
				for (size_t j = 0; j < sides2.size(); ++j) {
					match_and_identifiy(sides1[i], sides2[j], identifier);
				}
			}
		}
	}
}

template <class TElem>
void PeriodicBoundaryManager::match_and_identifiy(TElem* e1, TElem* e2,
		IIdentifier* i) {
	// try to use given identifier
	if (i != NULL)
		if (i->match(e1, e2)) {
			identifiy(e1, e2, i);
			return;
		}

	// fall back to member identifier
	if (!m_pIdentifier.get())
		UG_THROW("need an valid identifier.")

	if (m_pIdentifier->match(e1, e2)) {
		UG_LOG("match_and_ident: match found!!!\n")
		identifiy(e1, e2, m_pIdentifier.get());
	}
}

// is considered periodic if it is a master or a slave
template <class TElem>
bool PeriodicBoundaryManager::is_periodic(TElem* e) const {
	return get_periodic_status_accessor<TElem>()[e] != P_NOT_PERIODIC;
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
		if (master(*elem) == *elem) {
			Group<TElem>* g = group(*elem);
			UG_ASSERT(g, "group not valid")
			UG_LOG("group of " << (*elem)->reference_object_id() << "\tlevel: " <<
					m_pGrid->get_level(*elem) << "\tmaster: " <<
					GetGeometricObjectCenter(*m_pGrid, g->m_master) << "\tslaves: ");
			for (SlaveIter slave = g->get_slaves().begin();
					slave != g->get_slaves().end(); ++slave) {
				TElem* e = *slave;
				UG_LOG(GetGeometricObjectCenter(*m_pGrid, e) << ", ")
				UG_ASSERT(m_pGrid->get_level(*elem) == m_pGrid->get_level(e),
						"wrong level in group")
			}
			UG_LOG(std::endl)
		}
	}
}

template <class TElem, class TParent>
void PeriodicBoundaryManager::handle_creation(TElem* e, TParent* pParent,
		bool replacesParent) {

	typedef typename Group<TParent>::SlaveContainer ParentSlaveContainer;

	if (!pParent)
		return;

	if (!is_periodic(pParent))
		return;

	UG_ASSERT(m_pGrid,
			"PeriodicBoundaryManager has to operate on a grid. No grid assigned.");
	MultiGrid& mg = *m_pGrid;

	// if replacesParent == true, TParent == TElem
	if (replacesParent) {
		//	we only have to replace the parent entry.
		TElem* parent = dynamic_cast<TElem*>(pParent);
		if (parent) {
			if (is_periodic(parent)) {
				// if parent is master, set newly created item as new master
				if (is_master(parent)) {
					group(parent)->m_master = e;
					// set group of e
					set_group<TElem>(group(parent), e);
					// reset group of old parent
					set_group<TElem>(NULL, parent);
				} else { // slave
					// make newly created item a slave, since its parent is slave
					make_slave<TElem>(group(parent), e);
					set_group<TElem>(NULL, parent);
				}
			}
		} else
			UG_THROW("Can't replace an element with an element of another type.");
		// finished here
		return;
	}

	// attention: parentGroup may be null
	Group<TParent>* parentGroup = group<TParent>(pParent);

	if (is_master<TParent>(pParent)) {
		UG_LOG("parent is master, create group for e\n");
		// create new group for e, with e as master
		Group<TElem>* newGroup = new Group<TElem>(e);
		set_group(newGroup, e);

		// iterate over slaves of parent
		ParentSlaveContainer* parentSlaves = slaves(pParent);
		for (typename ParentSlaveContainer::iterator i_slave =
				parentSlaves->begin(); i_slave != parentSlaves->end();
				++i_slave) {
			TParent* parentSlave = *i_slave;

			// iterate over all children of parentSlave and make them slaves of e
			// only matching children of type TElem are considered.
			// note that not all children already have to exist at this point.
			// children which are created later on are handled by the 'else' section.
			for (size_t i_child = 0;
					i_child < mg.num_children<TElem>(parentSlave); ++i_child) {
				TElem* c = mg.get_child<TElem>(parentSlave, i_child);
				UG_ASSERT((!e->base_object_id() == VERTEX)
						  || (mg.num_children<TElem>(parentSlave) == 1),
						  "At most one 1 vertex-child is currently allowed");
				// We use a special case for vertices here, since the position
				// attachment has not yet been set (this is currently done after
				// all observers have been informed about the creation of the new
				// vertex). Note that this is no problem for edges or faces, since
				// associated vertices have been properly initialized at that time.
				if((e->base_object_id() == VERTEX) || (m_pIdentifier->match(c, e))) {
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
		if (parentGroup == NULL) {
			get_periodic_status_accessor<TElem>()[e] = P_SLAVE_MASTER_UNKNOWN;
			return;
		}

		TParent* parentMaster = parentGroup->m_master;
		bool master_found = false;
		for (size_t i_child = 0; i_child < mg.num_children<TElem>(parentMaster);
				++i_child) {
			TElem* c = mg.get_child<TElem>(parentMaster, i_child);
			UG_ASSERT((!e->base_object_id() == VERTEX)
					  || (mg.num_children<TElem>(parentMaster) == 1),
					  "At most one 1 vertex-child is currently allowed");
			// We use a special case for vertices here, since the position
			// attachment has not yet been set (this is currently done after
			// all observers have been informed about the creation of the new
			// vertex). Note that this is no problem for edges or faces, since
			// associated vertices have been properly initialized at that time.
			if((e->base_object_id() == VERTEX) || (m_pIdentifier->match(c, e))) {
				make_slave(group(c), e);
				master_found = true;
				break;
			}
		}

		// ATTENTION: If no master was found e currently can't be recognized as
		// a slave from other classes, as e.g. the DoF-Manager.
		// This, however, is important.
		// One should think about adding a master-/slave- flag and adjusting the
		// is_master(...) and is_slave(...) member methods. Furthermore it
		// should be documented, that a slave does not necessarily have a master
		// during grid adaption.
		// The add method in the dof-manager has to be adjusted, too. If a new slave
		// has been created, it should check whether its master already exists and
		// copy its dof-index from its master in this case. If the master doesn't
		// exist yet, the pseudo-index -1 should be assigned. The correct index
		// is then assigned later on during creation of the master, where a
		// dof-index is assigned to all existing slaves automatically.
		if (!master_found) {
			get_periodic_status_accessor<TElem>()[e] = P_SLAVE_MASTER_UNKNOWN;
		}
	}
}

/// handles deletion of element type
template <class TElem>
void PeriodicBoundaryManager::handle_deletion(TElem* e, TElem* replacedBy) {
//	UG_THROW("not tested.")
	if (!is_periodic(e))
		return;

	if (is_master(e)) {
		// if e is master, set its replacing element as new master
		if (replacedBy) {
			UG_ASSERT(!group(replacedBy),
					"replacing element is already in group")
			group(e)->m_master = replacedBy;
		} else {
			// todo whole group should be deleted
			remove_group(group(e));
		}
	} else { // slave
		Group<TElem>* g = group(e);
		bool removed = remove_slave(e);
		UG_ASSERT(removed, "slave not removed.")

		if (replacedBy)
			g->add_slave(replacedBy);
	}
}

template <class TElem>
void PeriodicBoundaryManager::make_slave(Group<TElem>* g, TElem* slave) {
	UG_ASSERT(m_pGrid->get_level(g->m_master) == m_pGrid->get_level(slave),
			"level of slave and group mismatch")
	g->add_slave(slave);
	set_group(g, slave);
}

template <class TElem>
bool PeriodicBoundaryManager::remove_slave(TElem* e) {
	if (slaves(e)) {
		typename Group<TElem>::SlaveContainer& s = *slaves(e);
		typename Group<TElem>::SlaveIterator pos = std::find(s.begin(), s.end(), e);
		if (pos != s.end()) {
			s.erase(pos);
			set_group<TElem>(NULL, e);
			// todo what if e was the last slave of of group?
			return true;
		}
	}
	return false;
}

template <class TElem>
void PeriodicBoundaryManager::remove_group(Group<TElem>* g) {
	if (!g)
		return;

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
		if (e && e != g0->m_master) {
			make_slave(g0, e);
		}
	}

	// make old master slave of group g1
	make_slave(g0, g1->m_master);

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
Grid::AttachmentAccessor<TElem,
		Attachment<PeriodicBoundaryManager::Group<TElem>*> >&
PeriodicBoundaryManager::get_group_accessor() {
	UG_THROW("not impled");
}

template <class TElem>
PeriodicBoundaryManager::Group<TElem>* PeriodicBoundaryManager::group(
		TElem* e) const {
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
		return;
	}

	if (g->m_master == e)
		get_periodic_status_accessor<TElem>()[e] = P_MASTER;
	else
		get_periodic_status_accessor<TElem>()[e] = P_SLAVE;
}

template <class TElem>
void PeriodicBoundaryManager::handle_creation_cast_wrapper(TElem* e,
		GeometricObject* pParent, bool replacesParent) {
	if (pParent) {
		switch (pParent->base_object_id()) {
		case VERTEX:
			handle_creation(e, static_cast<VertexBase*>(pParent),
					replacesParent);
			break;
		case EDGE:
			handle_creation(e, static_cast<EdgeBase*>(pParent), replacesParent);
			break;
		case FACE:
			handle_creation(e, static_cast<Face*>(pParent), replacesParent);
			break;
		// ignore volumes, as they are not meant to be periodic
		case VOLUME:
			break;
		default:
			UG_THROW("no handling for parent type: " << pParent->base_object_id())
		}
	} else {
		handle_creation(e, static_cast<VertexBase*>(NULL), replacesParent);
	}
}

#ifdef UG_DEBUG
/**
 * create all pairs of <master, slave> and insert them into a set to check for duplicates
 */
template <class TElem> void test(PeriodicBoundaryManager& PBM,
		GeometricObjectCollection& goc1, GeometricObjectCollection& goc2) {
	typedef typename ElementStorage<TElem>::SectionContainer::iterator gocIter;
	typedef typename std::pair<TElem*, TElem*> master_slave_pair;
	typedef typename std::set<master_slave_pair> MSSet;
	MSSet s;

	typedef typename PeriodicBoundaryManager::Group<TElem>::SlaveContainer SlaveContainer;
	typedef typename SlaveContainer::iterator SlaveIter;

	// subset 1
	for (gocIter iter = goc1.begin<TElem>(); iter != goc1.end<TElem>();
			++iter) {
		UG_ASSERT(PBM.is_periodic(*iter), "should be periodic now");

		if (PBM.master(*iter) == *iter) {
			SlaveContainer* slaves = PBM.slaves(*iter);
			UG_ASSERT(slaves, "master should have slaves")
			for (SlaveIter i = slaves->begin(); i != slaves->end(); ++i) {
				TElem* slave = *i;
				master_slave_pair p = std::make_pair(*iter, slave);
				bool inserted = (s.insert(p)).second;
				UG_ASSERT(inserted, "pair already exists");
			}
		}
	}
	// subset 1
	for (gocIter iter = goc2.begin<TElem>(); iter != goc2.end<TElem>();
			++iter) {
		UG_ASSERT(PBM.is_periodic(*iter), "should be periodic now");

		if (PBM.master(*iter) == *iter) {
			SlaveContainer* slaves = PBM.slaves(*iter);
			UG_ASSERT(slaves, "master should have slaves")
			for (SlaveIter i = slaves->begin(); i != slaves->end(); ++i) {
				TElem* slave = *i;
				master_slave_pair p = std::make_pair(*iter, slave);
				bool inserted = (s.insert(p)).second;
				UG_ASSERT(inserted, "pair already exists");
			}
		}
	}
}
#endif

template <class TDomain>
void IdentifySubsets(TDomain& dom, const char* sName1, const char* sName2) {
	// get subset handler from domain
	typedef typename TDomain::subset_handler_type subset_handler_type;

	subset_handler_type& sh = *dom.subset_handler();

	int si1 = sh.get_subset_index(sName1);
	int si2 = sh.get_subset_index(sName2);

	if (si1 == -1)
		UG_THROW("given subset name " << sName1 << " does not exist");
	if (si2 == -1)
		UG_THROW("given subset name " << sName2 << " does not exist");

	IdentifySubsets(dom, si1, si2);
}

/// performs geometric ident of periodic elements and master slave
template <class TDomain>
void IdentifySubsets(TDomain& dom, int sInd1, int sInd2) {
	if (sInd1 == -1 || sInd2 == -1) {
		UG_LOG("IdentifySubsets: at least one invalid subset given.\n")
		return;
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

	ParallelShiftIdentifier<position_accessor_type>* ident =
			new ParallelShiftIdentifier<position_accessor_type>(aaPos);
	pbm.set_identifier(ident);

	// shift vector between subsets
	position_type shift;

	// collect all geometric objects (even for all levels in case of multi grid)
	GeometricObjectCollection goc1 = sh.get_geometric_objects_in_subset(sInd1);
	GeometricObjectCollection goc2 = sh.get_geometric_objects_in_subset(sInd2);

	// map start type of recursion dependent to TDomain
	// in 3d start with faces, in 2d with edges, in 1d with vertices
	namespace mpl = boost::mpl;
	typedef		mpl::map<mpl::pair<Domain1d, VertexBase>,
						 mpl::pair<Domain2d, EdgeBase>,
						 mpl::pair<Domain3d, Face> > m;

	typedef typename mpl::at<m, TDomain>::type TElem;
	typedef typename ElementStorage<TElem>::SectionContainer::iterator gocIter;

	// calculate shift vector for top level
	position_type c1 = CalculateCenter(goc1.begin<TElem>(0), goc1.end<TElem>(0),
			aaPos);
	position_type c2 = CalculateCenter(goc2.begin<TElem>(0), goc2.end<TElem>(0),
			aaPos);

	VecSubtract(shift, c1, c2);
	ident->set_shift(shift);

	// for each level of multi grid. In case of simple Grid only one iteration
	for (size_t lvl = 0; lvl < goc1.num_levels(); lvl++) {
		// identify corresponding elements for second subset. A element is considered
		// to have symmetric element in second subset if there exists a shift vector between them.
		// todo use kd-tree for fast lookup of shifted elements
		for (gocIter iter1 = goc1.begin<TElem>(lvl);
				iter1 != goc1.end<TElem>(lvl); ++iter1)
			for (gocIter iter2 = goc2.begin<TElem>(lvl);
					iter2 != goc2.end<TElem>(lvl); ++iter2)
				pbm.match_and_identifiy(*iter1, *iter2, ident);
	}
#ifdef UG_DEBUG
	test<VertexBase>(pbm, goc1, goc2);
	test<EdgeBase>(pbm, goc1, goc2);
	test<Face>(pbm, goc1, goc2);
#endif
}
} // end namespace ug

#endif /* PERIODIC_IDENTIFIER_IMPL_HPP_ */

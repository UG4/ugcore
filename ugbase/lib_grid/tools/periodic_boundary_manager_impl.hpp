/*
 * periodic_identifier_impl.hpp
 *
 *  Created on: 26.11.2012
 *      Author: marscher
 */

#ifndef PERIODIC_IDENTIFIER_IMPL_HPP_
#define PERIODIC_IDENTIFIER_IMPL_HPP_

// include declarations
#include "periodic_boundary_manager.h"
#include "lib_grid/subset_handler.h"
#ifndef NDEBUG
	#include "lib_grid/algorithms/debug_util.h"
#else
	namespace ug {
	template <class TElem> TElem* GetGeometricObjectCenter(Grid& g, TElem* elem) {
		return elem;
	}}
#endif
#include "common/assert.h"

#include <boost/mpl/map.hpp>
#include <boost/mpl/at.hpp>

#include <algorithm>
#include <set>

namespace ug {

template<class TAAPos>
bool ParallelShiftIdentifier<TAAPos>::equals_shift(AttachmentType& diff)
{
	AttachmentType error;
	VecSubtract(error, diff, m_shift);

	if(std::abs(VecLength(error)) < 10E-8) return true;

	return false;
}

template<class TAAPos>
bool ParallelShiftIdentifier<TAAPos>::match(VertexBase* v1, VertexBase* v2)
{
	AttachmentType diff;
	VecSubtract(diff, m_aaPos[v1], m_aaPos[v2]);

	return equals_shift(diff);
}

template<class TAAPos>
bool ParallelShiftIdentifier<TAAPos>::match(EdgeBase* e1, EdgeBase* e2)
{
	AttachmentType c1 = CalculateCenter(e1, m_aaPos), c2 = CalculateCenter(e2,
			m_aaPos), diff;
	VecSubtract(diff, c1, c2);
	return equals_shift(diff);
}

template<class TAAPos>
bool ParallelShiftIdentifier<TAAPos>::match(Face* f1, Face* f2)
{
	AttachmentType c1 = CalculateCenter(f1, m_aaPos), c2 = CalculateCenter(f2,
			m_aaPos), diff;
	VecSubtract(diff, c1, c2);

	return equals_shift(diff);
}

template<class TAAPos>
bool ParallelShiftIdentifier<TAAPos>::match(Volume* v1, Volume* v2)
{
	AttachmentType c1 = CalculateCenter(v1, m_aaPos), c2 = CalculateCenter(v2,
			m_aaPos), diff;
	VecSubtract(diff, c1, c2);

	return equals_shift(diff);
}

template<class TElem>
void PeriodicBoundaryManager::identifiy(TElem* e1, TElem* e2,
		IIdentifier* identifier)
{
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
			} else if(m1 && m2) { // both are distinct masters
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
	if(identifier != NULL && TElem::HAS_SIDES) {
		container sides1, sides2;
		// collect sides and identify them
		m_pGrid->associated_elements<TElem>(sides1, e1);
		m_pGrid->associated_elements<TElem>(sides2, e2);
		UG_ASSERT(sides1.size() == sides2.size(),
				"sizes of periodic sub elements do not match!")
		for (size_t i = 0; i < sides1.size(); ++i) {
			for (size_t j = 0; j < sides2.size(); ++j) {
				// ensure elements are matching
				if(identifier->match(sides1[i], sides2[j])) {
					identifiy(sides1[i], sides2[j], identifier);
					break;
				}
			}
		}
	}
}

// is considered periodic if it is a master or a slave
template<class TElem>
bool PeriodicBoundaryManager::is_periodic(TElem* e) const
{
	return group(e) != NULL;
}

// gets master of e, may be null
template<class TElem>
TElem* PeriodicBoundaryManager::master(TElem* e) const
{
	if(group(e)) {
		return group(e)->m_master;
	}
	return NULL;
}

// gets slaves of e
template<class TElem>
std::list<TElem*>* PeriodicBoundaryManager::slaves(TElem* e) const
{
	if(group(e)) {
		return &group(e)->get_slaves();
	}
	return NULL;
}

template<class TElem>
void PeriodicBoundaryManager::print_identification() const
{
	typedef typename ElementStorage<TElem>::SectionContainer::iterator Iterator;
	typedef typename Group<TElem>::SlaveIterator SlaveIter;

	for (Iterator elem = m_pGrid->begin<TElem>(); elem != m_pGrid->end<TElem>();
			++elem) {
		// log masters and their slaves
		if(master(*elem) == *elem) {
			Group<TElem>* g = group(*elem);
			UG_ASSERT(g, "group not valid")
			UG_LOG("group: "<< g << " of " << (*elem)->reference_object_id()
					<<  "\t[master: " << GetGeometricObjectCenter(*m_pGrid, g->m_master) << "\tslaves: ");
			for (SlaveIter slave = g->get_slaves().begin();
					slave != g->get_slaves().end(); ++slave) {
				TElem* e = *slave;
				UG_LOG(GetGeometricObjectCenter(*m_pGrid, e) << ", ")
			}
			UG_LOG(std::endl)
		}
	}
}

template<class TElem, class TParent>
void PeriodicBoundaryManager::handle_creation(TElem* e,
		TParent* pParent, bool replacesParent) {
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
			UG_THROW(
					"Can't replace an element with an element of another type.");
		// finished here
		return;
	}

	if(pParent) {
		// if parent is not periodic, where are finished here
		if(!is_periodic<TParent>(pParent))
			return;

		// eg. elem => vertex, parent => vertex
		if(typeid(TElem) == typeid(TParent)) {
			TElem* parent = dynamic_cast<TElem*>(pParent);
			Group<TElem>* g = group(parent);
			UG_ASSERT(g, "group of parent not valid");
			make_slave(g, e);
		} else {
			// so parent is valid and periodic
			// get children of type TElem of parent
			uint children = m_pGrid->num_children<TElem>(pParent);
			UG_ASSERT(children > 0, "parent has no children.")
			bool child_found = false;
			for (uint i = 0; i < children; ++i) {
				TElem* child = m_pGrid->get_child<TElem>(pParent, i);
				UG_ASSERT(child, "child not valid");
				Group<TElem>* g = group<TElem>(child);
				if (g != NULL) {
					make_slave(g, e);
					child_found = true;
					break;
				}
			}
			// fixme this is assertion is always false!
			UG_ASSERT(child_found, "no periodic child found, which group could be copied.")
		}
	} else { // no parent
		// todo how to determine periodicity?
		UG_THROW("no parent, case not impled")
	}
}

/// handles deletion of element type
template<class TElem>
void PeriodicBoundaryManager::handle_deletion(TElem* e, TElem* replacedBy) {
	if(is_periodic(e)) {
		if(is_master(e)) {
			// if e is master, set its replacing element as new master
			if(replacedBy) {
				UG_ASSERT(!group(replacedBy), "replacing element is already in group")
				group(e)->m_master = replacedBy;
			} else {
				// todo whole group should be deleted
				remove_group(group(e));
			}
		} else { // slave
			Group<TElem>* g = group(e);
			bool removed = remove_slave(e);
			UG_ASSERT(removed, "slave not removed.")

			if(replacedBy)
				g->add_slave(replacedBy);
		}
	}
}

template <class TElem>
void PeriodicBoundaryManager::make_slave(Group<TElem>* g, TElem* slave)
{
	g->add_slave(slave);
	set_group(g, slave);
}

template <class TElem>
bool PeriodicBoundaryManager::remove_slave(TElem* e) {
	if(slaves(e)) {
		typename Group<TElem>::SlaveContainer& s = *slaves(e);
		typename Group<TElem>::SlaveIterator pos = std::find(s.begin(), s.end(), e);
		if(pos != s.end()) {
			s.erase(pos);
			set_group<TElem>(group(e), NULL);
			// todo what if e was the last slave of of group?
			return true;
		}
	}
	return false;
}

template <class TElem>
void PeriodicBoundaryManager::remove_group(Group<TElem>* g) {
	if(!g)
		return;

	// reset group pointers of all group members to NULL
	set_group<TElem>(NULL, g->m_master);
	typename Group<TElem>::SlaveContainer& s = g->get_slaves();
	typename Group<TElem>::SlaveIterator iter;
	for(iter = s.begin(); iter != s.end(); ++iter) {
		set_group<TElem>(NULL, *iter);
	}

	delete g;
}


/**
 * merges g1 in g0 and deletes g1 afterwards
 */
template<class TElem>
void PeriodicBoundaryManager::merge_groups(Group<TElem>* g0,
		Group<TElem>* g1)
{
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
	if(master(e) == NULL) {
		return false;
	}

	return master(e) != e;
}

template <class TElem>
bool PeriodicBoundaryManager::is_master(TElem* e) const {
	if(master(e) == e) {
		return true;
	}

	return false;
}

template <class TElem>
const Grid::AttachmentAccessor<TElem, Attachment<PeriodicBoundaryManager::Group<TElem>* > >&
PeriodicBoundaryManager::get_group_accessor() const {
	UG_THROW("not impled");
}

template<class TElem>
PeriodicBoundaryManager::Group<TElem>* PeriodicBoundaryManager::group(
		TElem* e) const
{
	return get_group_accessor<TElem>()[e];
}

template<class TElem>
void PeriodicBoundaryManager::set_group(Group<TElem>* g, TElem* e)
{
	const_cast<Grid::AttachmentAccessor<TElem, Attachment<Group<TElem>* > >&>
		(get_group_accessor<TElem>())[e] = g;
}

template <class TElem>
void PeriodicBoundaryManager::handle_creation_cast_wrapper(TElem* e, GeometricObject* pParent, bool replacesParent) {
	if(pParent) {
		switch(pParent->base_object_id()) {
		case VERTEX: handle_creation(e, static_cast<VertexBase*>(pParent), replacesParent); break;
		case EDGE: handle_creation(e, static_cast<EdgeBase*>(pParent), replacesParent); break;
		case FACE: handle_creation(e, static_cast<Face*>(pParent), replacesParent); break;
		default:
			UG_THROW("no handling for parent type: " << pParent->base_object_id())
		}
	} else {
		handle_creation(e, static_cast<VertexBase*>(NULL), replacesParent);
	}
}

#ifndef NDEBUG
/**
 * create all pairs of <master, slave> and insert them into a set to check for duplicates
 */
template<class TElem> void test(PeriodicBoundaryManager& PBM,
		GeometricObjectCollection& goc1, GeometricObjectCollection& goc2)
{
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

		if(PBM.master(*iter) == *iter) {
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

		if(PBM.master(*iter) == *iter) {
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

	if(si1 == -1)
		UG_THROW("given subset name " << sName1 << " does not exist");
	if(si2 == -1)
		UG_THROW("given subset name " << sName2 << " does not exist");

	IdentifySubsets(dom, si1, si2);
}

/// performs geometric ident of periodic elements and master slave
template<class TDomain>
void IdentifySubsets(TDomain& dom, int sInd1, int sInd2)
{
	if(sInd1 == -1 || sInd2 == -1) {
		UG_LOG("IdentifySubsets: at least one invalid subset given.\n")
		return;
	}

	// ensure grid has support for periodic boundaries
	if(!dom.grid()->has_periodic_boundaries())
	{
		dom.grid()->set_periodic_boundaries(true);
	}

	PeriodicBoundaryManager& PI = *dom.grid()->periodic_boundary_manager();

//	UG_ASSERT(!dom.is_adaptive(), "adaptive domains currently not supported!");

	typedef typename TDomain::position_type position_type;
	typedef typename TDomain::position_accessor_type position_accessor_type;

	// get subset handler from domain
	typedef typename TDomain::subset_handler_type subset_handler_type;

	subset_handler_type& sh = *dom.subset_handler();

	// get aaPos from domain
	position_accessor_type& aaPos = dom.position_accessor();

	ParallelShiftIdentifier<position_accessor_type> ident(aaPos);

	// some vectors
	position_type c1, c2, shift, diff, error;

	// collect all geometric objects (even for all levels in case of multi grid)
	GeometricObjectCollection goc1 = sh.get_geometric_objects_in_subset(sInd1);
	GeometricObjectCollection goc2 = sh.get_geometric_objects_in_subset(sInd2);

	// map start type of recursion dependent to TDomain
	// in 3d start with faces, in 2d with edges, in 1d with vertices
	namespace mpl = boost::mpl;
	typedef mpl::map<mpl::pair<Domain1d, VertexBase>,
					  mpl::pair<Domain2d, EdgeBase>,
					  mpl::pair<Domain3d, Face> > m;

	typedef typename mpl::at<m, TDomain>::type TElem;
	typedef typename ElementStorage<TElem>::SectionContainer::iterator gocIter;

	// for each level of multi grid. In case of simple Grid only one iteration
	for (size_t lvl = 0; lvl < goc1.num_levels(); lvl++) {

		// calculate shift vector for current level
		c1 = CalculateCenter(goc1.begin<TElem>(lvl), goc1.end<TElem>(lvl),
				aaPos);
		c2 = CalculateCenter(goc2.begin<TElem>(lvl), goc2.end<TElem>(lvl),
				aaPos);

		// calculate shift vector for objects
		VecSubtract(shift, c1, c2);
		ident.set_shift(shift);

		// identify corresponding elements for second subset. A element is considered
		// to have symmetric element in second subset if there exists a shift vector between them.
		// todo use kd-tree for fast lookup of shifted elements
		for (gocIter iter1 = goc1.begin<TElem>(lvl);
				iter1 != goc1.end<TElem>(lvl); ++iter1) {
			for (gocIter iter2 = goc2.begin<TElem>(lvl);
					iter2 != goc2.end<TElem>(lvl); ++iter2) {
				if(ident.match(*iter1, *iter2)) PI.identifiy(*iter1, *iter2,
						&ident);
			}
		}
	}
#ifndef NDEBUG
	test<VertexBase>(PI, goc1, goc2);
	test<EdgeBase>(PI, goc1, goc2);
	test<Face>(PI, goc1, goc2);
#endif
}
} // end namespace ug

#endif /* PERIODIC_IDENTIFIER_IMPL_HPP_ */

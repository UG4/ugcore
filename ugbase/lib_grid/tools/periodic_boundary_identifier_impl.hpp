/*
 * periodic_identifier_impl.hpp
 *
 *  Created on: 26.11.2012
 *      Author: marscher
 */

#ifndef PERIODIC_IDENTIFIER_IMPL_HPP_
#define PERIODIC_IDENTIFIER_IMPL_HPP_

// include declarations
#include "periodic_boundary_identifier.h"
#include "lib_grid/subset_handler.h"
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

// initial call with TElem = ug::Volume to assure all associated sub elements are identified
template<class TElem>
void PeriodicBoundaryIdentifier::identifiy(TElem* e1, TElem* e2,
		IIdentifier* identifier)
{
	UG_LOG("ident\n")
	typedef typename Grid::traits<typename TElem::side>::secure_container container;

	// determine master
	TElem* m1 = master(e1);
	TElem* m2 = master(e2);

	// assign groups
	if(m1 && !m2) // m1 is master, so e2 will be its slave
	{
		UG_LOG("e1 is master, e2 not\n")
		slaves(m1)->push_back(e2);
	}
	else if(m2 && !m1) // m2 is master, so e1 will be its slave
	{
		UG_LOG("e2 is master, e1 not\n")
		slaves(m2)->push_back(e1);
	}
	else if(m1 && m2) // both are masters, merge them.
	{
		UG_LOG("both are masters\n")
		merge_groups(group(e1), group(e2));
	}
	else // nobody is master
	{
		UG_LOG("nobody is master\n")
		// create new group with e1 as master and e2 as slave
		Group<TElem>* g = new Group<TElem>(e1);
		g->add_slave(e2);

		set_group(g, e1);
		set_group(g, e2);
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
bool PeriodicBoundaryIdentifier::is_periodic(TElem* e) const
{
	return group(e) != NULL;
}

// gets master of e, may be null
template<class TElem>
TElem* PeriodicBoundaryIdentifier::master(TElem* e) const
{
	if(group(e)) {
		return group(e)->m_master;
	}
	return NULL;
}

// gets slaves of e
template<class TElem>
std::list<TElem*>* PeriodicBoundaryIdentifier::slaves(TElem* e) const
{
	if(group(e)) {
		return &group(e)->get_slaves();
	}
	return NULL;
}

template<class TElem>
void PeriodicBoundaryIdentifier::print_identification() const
{
	typedef typename ElementStorage<TElem>::SectionContainer::iterator Iterator;
	typedef typename Group<TElem>::SlaveIterator SlaveIter;
	for (Iterator iter = m_pGrid->begin<TElem>(); iter != m_pGrid->end<TElem>();
			++iter) {
		// log masters and their slaves
		if(master(*iter)) {
			Group<TElem>* g = group(*iter);
			UG_LOG("group: [master: " << g->m_master << "\tslaves: ");
			for (SlaveIter i = g->get_slaves().begin();
					i != g->get_slaves().end(); ++i) {
				UG_LOG(*i << ", ")
			}
			UG_LOG(std::endl)
		}
	}
}

/**
 * merges g1 in g0 and deletes g1 afterwards
 */
template<class TElem>
void PeriodicBoundaryIdentifier::merge_groups(Group<TElem>* g0,
		Group<TElem>* g1)
{
	UG_LOG("merge called" << std::endl)
	typedef typename Group<TElem>::SlaveContainer SlaveContainer;
	typedef typename Group<TElem>::SlaveIterator SlaveIterator;

	SlaveContainer& slaves_g0 = g0->get_slaves();
	SlaveContainer& slaves_g1 = g1->get_slaves();

	// insert slaves of g1 at the end of slaves of g0
	slaves_g0.insert(slaves_g0.end(), slaves_g1.begin(), slaves_g1.end());
	// insert prior master of g1 to slaves of g0
	slaves_g0.push_back(g1->m_master);

	// set group pointer to g0 for all slaves of g1
	for (SlaveIterator iter = slaves_g1.begin(); iter != slaves_g1.end(); ++iter) {
		set_group(g0, *iter);
	}
	// remove old group
	delete g1;
}

template<class TElem>
PeriodicBoundaryIdentifier::Group<TElem>* PeriodicBoundaryIdentifier::group(
		TElem* e) const
{
	UG_THROW("not impled");
}

template<>
PeriodicBoundaryIdentifier::Group<VertexBase>* PeriodicBoundaryIdentifier::group(
		VertexBase* e) const
{
	return m_aaGroupVRT[e];
}

template<>
PeriodicBoundaryIdentifier::Group<EdgeBase>* PeriodicBoundaryIdentifier::group(
		EdgeBase* e) const
{
	return m_aaGroupEDG[e];
}

template<>
PeriodicBoundaryIdentifier::Group<Face>* PeriodicBoundaryIdentifier::group(
		Face* e) const
{
	return m_aaGroupFCE[e];
}

template<>
PeriodicBoundaryIdentifier::Group<Volume>* PeriodicBoundaryIdentifier::group(
		Volume* e) const
{
	return m_aaGroupVOL[e];
}

template<class TElem>
void PeriodicBoundaryIdentifier::set_group(Group<TElem>* g, TElem* e)
{
	UG_THROW("not impled")
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<VertexBase>* g, VertexBase* v)
{
	m_aaGroupVRT[v] = g;
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<EdgeBase>* g, EdgeBase* e)
{
	m_aaGroupEDG[e] = g;
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<Face>* g, Face* f)
{
	m_aaGroupFCE[f] = g;
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<Volume>* g, Volume* v)
{
	m_aaGroupVOL[v] = g;
}

void PeriodicBoundaryIdentifier::set_grid(Grid* g)
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

PeriodicBoundaryIdentifier::~PeriodicBoundaryIdentifier()
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
	// fixme no groups for vols in 3d?
	for (VolumeIterator iter = m_pGrid->begin<Volume>();
			iter != m_pGrid->end<Volume>(); ++iter) {
		if(master(*iter)) delete m_aaGroupVOL[*iter];
	}

	// set grid to NULL to detach groups from grid
	set_grid (NULL);
}

template<class TElem> void test(PeriodicBoundaryIdentifier& PI,
		GeometricObjectCollection& goc1, GeometricObjectCollection& goc2)
{
	UG_LOG("test\n")
	typedef typename ElementStorage<TElem>::SectionContainer::iterator gocIter;
	typedef typename std::pair<TElem*, TElem*> master_slave_pair;
	typedef typename std::set<master_slave_pair> MSSet;
	MSSet s;

	typedef typename PeriodicBoundaryIdentifier::Group<TElem>::SlaveContainer SlaveContainer;
	typedef typename SlaveContainer::iterator SlaveIter;
	// subset 1
	for (gocIter iter = goc1.begin<TElem>(); iter != goc1.end<TElem>();
			++iter) {
		UG_ASSERT(PI.is_periodic(*iter), "should be periodic now");

		if(PI.master(*iter) == *iter) {
			SlaveContainer* slaves = PI.slaves(*iter);
			UG_ASSERT(slaves, "master should have slaves")
			for (SlaveIter i = slaves->begin(); i != slaves->end(); ++i) {
				TElem* slave = *i;
				master_slave_pair p = std::make_pair(*iter, slave);
				bool inserted = (s.insert(p)).second;
				UG_ASSERT(inserted, "pair already exists");
			}
		}
	}

//	// subset 2
//	for (gocIter iter = goc2.begin<TElem>(); iter != goc2.end<TElem>();
//			++iter) {
//		if(PI.master(*iter) == *iter) {
//			std::vector<TElem*>* slaves = PI.slaves(*iter);
//			UG_ASSERT(slaves, "master should have slaves")
//			for (uint i = 0; i < slaves->size(); i++) {
//				TElem* slave = slaves->at(i);
//				master_slave_pair p = std::make_pair(*iter, slave);
//				bool inserted = (s.insert(p)).second;
//				UG_ASSERT(inserted, "pair already exists");
//			}
//		}
//	}

}

// performs geometric ident of periodic elements and master slave
template<class TDomain>
void IdentifySubsets(TDomain& dom, PeriodicBoundaryIdentifier& PI, int sInd1,
		int sInd2)
{

//	UG_ASSERT(!dom.is_adaptive(), "adaptive domains currently not supported!");

	typedef typename TDomain::position_type position_type;
	typedef typename TDomain::position_accessor_type position_accessor_type;

	// get subset handler from domain
	typedef typename TDomain::subset_handler_type subset_handler_type;

	subset_handler_type& sh = *dom.subset_handler();

	// get aaPos from domain
	position_accessor_type& aaPos = dom.position_accessor();

	ParallelShiftIdentifier<position_accessor_type> ident(aaPos);

	// set grid of periodic identifier
	PI.set_grid(sh.grid());

	// some vectors
	position_type c1, c2, shift, diff, error;

	// collect all geometric objects (even for all levels in case of multi grid)
	GeometricObjectCollection goc1 = sh.get_geometric_objects_in_subset(sInd1);
	GeometricObjectCollection goc2 = sh.get_geometric_objects_in_subset(sInd2);

	// map start type of recursion dependant to TDomain
	namespace mpl = boost::mpl;
	typedef mpl::map<mpl::pair<Domain1d, VertexBase>,
			mpl::pair<Domain2d, EdgeBase>, mpl::pair<Domain3d, Face> > m;

	// in 3d start with faces
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
				// fixme this leads to duplicate ident on this element level
				if(ident.match(*iter1, *iter2)) PI.identifiy(*iter1, *iter2,
						&ident);
			}
		}
	}
//#ifdef DEBUG
	test<VertexBase>(PI, goc1, goc2);
	test<EdgeBase>(PI, goc1, goc2);
	test<Face>(PI, goc1, goc2);
//#endif
}
} // end namespace ug

#endif /* PERIODIC_IDENTIFIER_IMPL_HPP_ */

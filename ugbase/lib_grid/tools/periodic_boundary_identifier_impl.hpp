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

#include <algorithm>

namespace ug {

bool PeriodicBoundaryIdentifier::IIdentifier::match(VertexBase* v1, VertexBase* v2) {
	// todo impl
}

// initial call with TElem = ug::Volume to assure all associated sub elements are identified
template<class TElem>
void PeriodicBoundaryIdentifier::identifiy(TElem* e1, TElem* e2, IIdentifier* identifier)
{

	typedef typename Grid::traits<typename TElem::side>::secure_container container;

	// determine master
	TElem* m1 = master(e1);
	TElem* m2 = master(e2);

	// assign groups
	if(m1 && !m2) // m1 is master, so e2 will be its slave
	{
		slaves(m1)->push_back(e2);
	}
	else if(m2 && !m1) // m2 is master, so e1 will be its slave
	{
		slaves(m2)->push_back(e1);
	}
	else if(m1 && m2) // both are masters, merge them.
	{
		Group<TElem>* g1 = group(e1);
		Group<TElem>* g2 = group(e2);
		merge_groups(g1, g2);
	}
	else // nobody is master
	{
		// create new group with e1 as master and e2 as slave
		Group<TElem>* g = new Group<TElem>(e1);
		g->add_slave(e2);

		set_group(g, e1);
		set_group(g, e2);
	}

	// while elements have sides recursively identify sub type elements
	if(TElem::HAS_SIDES) {
		container sides1, sides2;
		// collect sides and identify them
		m_pGrid->associated_elements<TElem>(sides1, e1);
		m_pGrid->associated_elements<TElem>(sides2, e2);
		UG_ASSERT(sides1.size() == sides2.size(),
				"sizes of periodic sub elements do not match!")
		for (size_t i = 0; i < sides1.size(); ++i) {
			// fixme use iidentification match stuff
			identifiy(sides1[i], sides2[i]);
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
std::vector<TElem*>* PeriodicBoundaryIdentifier::slaves(TElem* e) const
{
	if(group(e)) {
		return &group(e)->get_slaves();
	}
	return NULL;
}

/**
 * merges g1 in g0 and deletes g1 afterwards
 */
template<class TElem>
void PeriodicBoundaryIdentifier::merge_groups(Group<TElem>* g0,
		Group<TElem>* g1)
{
	std::vector<TElem*>& slaves_g0 = g0->get_slaves(), slaves_g1 =
			g1->get_slaves();
	typedef typename std::vector<TElem*>::iterator iterator;

	// insert slaves of g1 at the end of slaves of g0
	slaves_g0.insert(slaves_g0.end(), slaves_g1.begin(), slaves_g1.end());
	// insert prior master of g1 to slaves of g0
	g0->add_slave(g1->m_master);

	// set Group for all elements of g1 to g0
	for (iterator iter = slaves_g1.begin(); iter != slaves_g1.end(); ++iter) {
		set_group(g0, *iter);
	}

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
	return m_aaGroupInfoVRT[e];
}

template<>
PeriodicBoundaryIdentifier::Group<EdgeBase>* PeriodicBoundaryIdentifier::group(
		EdgeBase* e) const
{
	return m_aaGroupInfoEDG[e];
}

template<>
PeriodicBoundaryIdentifier::Group<Face>* PeriodicBoundaryIdentifier::group(
		Face* e) const
{
	return m_aaGroupInfoFCE[e];
}

template<>
PeriodicBoundaryIdentifier::Group<Volume>* PeriodicBoundaryIdentifier::group(
		Volume* e) const
{
	return m_aaGroupInfoVOL[e];
}

template<class TElem>
void PeriodicBoundaryIdentifier::set_group(Group<TElem>* g, TElem* e)
{
	UG_THROW("not impled")
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<VertexBase>* g, VertexBase* v)
{
	m_aaGroupInfoVRT[v] = g;
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<EdgeBase>* g, EdgeBase* e)
{
	m_aaGroupInfoEDG[e] = g;
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<Face>* g, Face* f)
{
	m_aaGroupInfoFCE[f] = g;
}

template<>
void PeriodicBoundaryIdentifier::set_group(Group<Volume>* g, Volume* vol)
{
	m_aaGroupInfoVOL[vol] = g;
}

void PeriodicBoundaryIdentifier::set_grid(Grid* g)
{
	if(g != NULL) m_pGrid = g;

	// todo avoid creation here
	// group attachments
	Attachment<Group<VertexBase>* > aGroupVRT;
	Attachment<Group<EdgeBase>* > aGroupEDG;
	Attachment<Group<Face>* > aGroupFCE;
	Attachment<Group<Volume>* > aGroupVOL;

	// if grid changes, attach groups
	if(m_pGrid != NULL && m_pGrid != g) {
		// access grid with those attachments
		m_aaGroupInfoVRT.access(*m_pGrid, aGroupVRT);
		m_aaGroupInfoEDG.access(*m_pGrid, aGroupEDG);
		m_aaGroupInfoFCE.access(*m_pGrid, aGroupFCE);
		m_aaGroupInfoVOL.access(*m_pGrid, aGroupVOL);
	}

	// detach groups
	if(g == NULL)
	{
		m_pGrid->detach_from_vertices(aGroupVRT);
		m_pGrid->detach_from_edges(aGroupEDG);
		m_pGrid->detach_from_faces(aGroupFCE);
		m_pGrid->detach_from_volumes(aGroupVOL);
	}
}

PeriodicBoundaryIdentifier::~PeriodicBoundaryIdentifier() {
	// delete group instances of all masters for all elements
	for (VertexBaseIterator iter = m_pGrid->begin<VertexBase>();
			iter != m_pGrid->end<VertexBase>(); ++iter) {
		if(master(*iter))
			delete m_aaGroupInfoVRT[*iter];
	}

	for (EdgeBaseIterator iter = m_pGrid->begin<EdgeBase>();
			iter != m_pGrid->end <EdgeBase>(); ++iter) {
		if(master(*iter))
			delete m_aaGroupInfoEDG[*iter];
	}

	for (FaceIterator iter = m_pGrid->begin<Face>();
			iter != m_pGrid->end <Face>(); ++iter) {
		if(master(*iter))
			delete m_aaGroupInfoFCE[*iter];
	}

	for (VolumeIterator iter = m_pGrid->begin<Volume>();
			iter != m_pGrid->end <Volume>(); ++iter) {
		if(master(*iter))
			delete m_aaGroupInfoVOL[*iter];
	}

	// set grid to NULL to detach groups from grid
	set_grid(NULL);
}

// performs geometric ident of periodic elements and master slave
template<class TDomain>
void IdentifySubsets(TDomain& dom, PeriodicBoundaryIdentifier& PI, int sInd1,
		int sInd2)
{

	UG_ASSERT(!dom.is_adaptive(), "adaptive domains currently not supported!");

	typedef typename TDomain::position_type position_type;
	typedef typename TDomain::position_accessor_type position_accessor_type;

	// get subset handler from domain
	typedef typename TDomain::subset_handler_type subset_handler_type;

	subset_handler_type& sh = *dom.subset_handler();

	// get aaPos from domain
	position_accessor_type& aaPos = dom.position_accessor();

	// set grid of periodic identifier
	PI.set_grid(sh.grid());

	// some vectors
	position_type c1, c2, shift, diff, error;

	// collect all geometric objects (even for all levels in case of multi grid)
	GeometricObjectCollection goc1 = sh.get_geometric_objects_in_subset(sInd1);
	GeometricObjectCollection goc2 = sh.get_geometric_objects_in_subset(sInd2);

	// start identification of volumes, which will be recursively identified to its vertices
	typedef Volume TElem;
	// for each level of multi grid. In case of simple Grid only one iteration
	for (size_t lvl = 0; lvl < goc1.num_levels(); lvl++) {

		c1 = CalculateCenter(goc1.begin<TElem>(lvl), goc1.end<TElem>(lvl),
				aaPos);
		c2 = CalculateCenter(goc2.begin<TElem>(lvl), goc2.end<TElem>(lvl),
				aaPos);

		// calculate shift vector for objects
		VecSubtract(shift, c1, c2);

		// identify corresponding elements for second subset. A element is considered
		// to have symmetric element in second subset if there exists a shift vector between them.
		// todo use kd-tree for fast lookup of shifted elements
		typedef ElementStorage<TElem>::SectionContainer::iterator gocIter;
		for (gocIter iter1 = goc1.begin<TElem>(lvl);
				iter1 != goc1.end<TElem>(lvl); ++iter1) {
			c1 = CalculateCenter(*iter1, aaPos);
			for (gocIter iter2 = goc2.begin<TElem>(lvl);
					iter2 != goc2.end<TElem>(lvl); ++iter2) {
				c2 = CalculateCenter(*iter2, aaPos);
				VecSubtract(diff, c1, c2);
				VecSubtract(error, shift, diff);

				if(std::abs(VecLength(error)) < 10E-8) {
					PI.identifiy(*iter1, *iter2);
					break;
				}
			}
		}
	}
}

} // end namespace ug

#endif /* PERIODIC_IDENTIFIER_IMPL_HPP_ */

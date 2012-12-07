/*
 * periodic_boundary_manager.h
 *
 *  Created on: 26.11.2012
 *      Author: marscher
 */

#ifndef PERIODIC_IDENTIFIER_H_
#define PERIODIC_IDENTIFIER_H_

#include "lib_grid/grid/grid.h"
#include "lib_grid/grid/geometric_base_objects.h"
#include "lib_disc/domain.h"

#include <vector>
#include <iostream>

namespace ug {

/// Interface to match periodic geometric elements
/**
 *
 */
class IIdentifier {
public:
	virtual bool match(VertexBase*, VertexBase*) = 0;
	virtual bool match(EdgeBase*, EdgeBase*) = 0;
	virtual bool match(Face*, Face*) = 0;
	virtual bool match(Volume*, Volume*) = 0;
};

/// This class matches geometric elements which are parallel translated.
/**
 * Usage: class needs to be instantiated with the position attachment used on the Domain.
 * Before using any match methods, the translation vector needs to be set with set_shift()
 * \tparam <TPosAA>{class needs to be instantiated with the
 * position attachment used on the Domain.}
 */
template<class TPosAA> class ParallelShiftIdentifier: public IIdentifier {
public:
	virtual bool match(VertexBase*, VertexBase*);
	virtual bool match(EdgeBase*, EdgeBase*);
	virtual bool match(Face*, Face*);
	virtual bool match(Volume*, Volume*);

	typedef typename TPosAA::ValueType AttachmentType;
	ParallelShiftIdentifier(TPosAA& aa) : m_aaPos(aa) {}
	void set_shift(AttachmentType& shift) {m_shift = shift;}
protected:
	AttachmentType m_shift;
	TPosAA& m_aaPos;
	bool equals_shift(AttachmentType& diff);
};

///
/**
 *
 */
class PeriodicBoundaryManager {
public:
	PeriodicBoundaryManager() : m_pGrid(NULL) {}
	~PeriodicBoundaryManager();

	// sets grid and inits group info attachment accessors
	void set_grid(Grid* g);

	///
	/**
	 *
	 */
	template<class TElem> void identifiy(TElem* e1, TElem* e2, IIdentifier* i = NULL);
	template<class TElem> bool is_periodic(TElem* e) const;
	template<class TElem> bool is_slave(TElem*) const;
	template<class TElem> bool is_master(TElem*) const;
	template<class TElem> TElem* master(TElem* e) const;
	// todo use Group<TElem>::SlaveContainer as return type
	template<class TElem> std::list<TElem*>* slaves(TElem* e) const;
	template<class TElem> void print_identification() const;

	/**
	 * A Group instance holds a master of type TElem and several children.
	 */
	template <class TElem> class Group {
	public:
		typedef typename std::list<TElem*> SlaveContainer;
		typedef typename SlaveContainer::iterator SlaveIterator;

		Group(TElem* m = NULL) : m_master(m) {}

		void add_slave(TElem* e) {
			UG_ASSERT(e, "add_slave: slave not valid");
			UG_ASSERT(e != m_master, "element already master!");
			m_slaves.push_back(e); }

		SlaveContainer& get_slaves() { return m_slaves; }
		TElem* m_master;

	protected:
		SlaveContainer m_slaves;
	};

protected:
	// no copy construction allowed
	PeriodicBoundaryManager(const PeriodicBoundaryManager&)	{}

	// grid instance we operate on
	Grid* m_pGrid;

	/// attachment accessors for Groups
	Grid::AttachmentAccessor<VertexBase, Attachment<Group<VertexBase>* > > m_aaGroupVRT;
	Grid::AttachmentAccessor<EdgeBase, Attachment<Group<EdgeBase>* > > m_aaGroupEDG;
	Grid::AttachmentAccessor<Face, Attachment<Group<Face>* > > m_aaGroupFCE;
	Grid::AttachmentAccessor<Volume, Attachment<Group<Volume>* > > m_aaGroupVOL;

	template <class TElem> void make_slave(Group<TElem>* g, TElem* slave);
	template <class TElem> void merge_groups(Group<TElem>* g0, Group<TElem>* g1);

	// get typed attachment accessor for group attachment
	template <class TElem>
	const Grid::AttachmentAccessor<TElem, Attachment<Group<TElem>* > >&
	get_group_accessor() const;

	// gets group of given element e
	template <class TElem> Group<TElem>* group(TElem* e) const;
	// set group attachment to element e
	template <class TElem> void set_group(Group<TElem>* g, TElem* e);
};

/**
 * \brief identifies subset 1 with subset 2. If the grid of given domain has no
 * periodic boundary manager attached, one will be created.
 *
 * \param dom Domain the periodic boundary should be defined on
 * \param sInd1 subset index which elements should be identified with elements from
 * those of sInd2
 * \param sInd2 \see{sInd1}
 */
template <class TDomain>
void IdentifySubsets(TDomain& dom, int sInd1, int sInd2);

/**
 * \brief identifies subset 1 with subset 2. If the grid of given domain has no
 * periodic boundary manager attached, one will be created.
 *
 * \param dom Domain the periodic boundary should be defined on
 * \param sName1 subset name which elements should be identified with elements from
 * those of sName2
 * \param sName2 \see {sName1}
 */
template <class TDomain>
void IdentifySubsets(TDomain& dom, const char* sName1, const char* sName2);
} // end of namespace ug

// include implementation
#include "periodic_boundary_manager_impl.hpp"

#endif /* PERIODIC_IDENTIFIER_H_ */

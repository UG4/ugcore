/*
 * periodic_identifier.h
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

class IIdentifier {
public:
	virtual bool match(VertexBase*, VertexBase*){UG_THROW("not impled");}
	virtual bool match(EdgeBase*, EdgeBase*){UG_THROW("not impled");}
	virtual bool match(Face*, Face*){UG_THROW("not impled");}
	virtual bool match(Volume*, Volume*){UG_THROW("not impled");}
};

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


/**
 *
 */
class PeriodicBoundaryIdentifier {
public:
	PeriodicBoundaryIdentifier() : m_pGrid(NULL) {}
	~PeriodicBoundaryIdentifier();

	// sets grid and inits group info attachment accessors
	void set_grid(Grid* g);

	template<class TElem> void identifiy(TElem* e1, TElem* e2, IIdentifier* i = NULL);
	template<class TElem> bool is_periodic(TElem* e) const;
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
		~Group() {UG_LOG("group of master " << m_master << " destroyed\n")}

		void add_slave(TElem* e) { UG_ASSERT(e != m_master, "element already master!"); m_slaves.push_back(e); }
		SlaveContainer& get_slaves() { return m_slaves; }
		TElem* m_master;

	protected:
		SlaveContainer m_slaves;
	};

protected:
	// grid instance we operate on
	Grid* m_pGrid;

	/// attachment accessors for Groups
	Grid::AttachmentAccessor<VertexBase, Attachment<Group<VertexBase>* > > m_aaGroupVRT;
	Grid::AttachmentAccessor<EdgeBase, Attachment<Group<EdgeBase>* > > m_aaGroupEDG;
	Grid::AttachmentAccessor<Face, Attachment<Group<Face>* > > m_aaGroupFCE;
	Grid::AttachmentAccessor<Volume, Attachment<Group<Volume>* > > m_aaGroupVOL;

	template <class TElem> void merge_groups(Group<TElem>* g0, Group<TElem>* g1);

	// gets group of given element
	template <class TElem> Group<TElem>* group(TElem* e) const;

	template <class TElem> void set_group(Group<TElem>* g, TElem* e);
};

/**
 * \brief
 *
 * \param dom
 * \param PI
 * \param sInd1
 * \param sInd2
 */
template <class TDomain>
void IdentifySubsets(TDomain& dom, PeriodicBoundaryIdentifier& PI, int sInd1, int sInd2);

} // end of namespace ug

// include implementation
#include "periodic_boundary_identifier_impl.hpp"

#endif /* PERIODIC_IDENTIFIER_H_ */

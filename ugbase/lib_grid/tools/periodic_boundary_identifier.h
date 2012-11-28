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

namespace ug {

/**
 *
 */
class PeriodicBoundaryIdentifier {

public:
//	PeriodicBoundaryIdentifier() : m_pGrid(NULL) {}
	~PeriodicBoundaryIdentifier();

	// sets grid and inits group info attachment accessors
	void set_grid(Grid* g);

	class IIdentifier {
	public:
		virtual bool match(VertexBase*, VertexBase*);
		virtual bool match(EdgeBase*, EdgeBase*);
		virtual bool match(Face*, Face*);
		virtual bool match(Volume*, Volume*);
	};

	template<class TElem> void identifiy(TElem* e1, TElem* e2, IIdentifier* i = NULL);
	template<class TElem> bool is_periodic(TElem* e) const;
	template<class TElem> TElem* master(TElem* e) const;
	template<class TElem> std::vector<TElem*>* slaves(TElem* e) const;

protected:
	// grid instance we operate on, set by set_
	Grid* m_pGrid;

	/**
	 * A Group instance holds a master of type TElem and several children.
	 */
	template <class TElem> class Group {
	public:
		Group(TElem* m = NULL) : m_master(m) {}
		void add_slave(TElem* e) { UG_ASSERT(e != m_master, "duplicate master!"); m_slaves.push_back(e); }
		std::vector<TElem*>& get_slaves() { return m_slaves; }
		TElem* m_master;
	protected:
		std::vector<TElem*> m_slaves;
	};

	Grid::AttachmentAccessor<VertexBase, Attachment<Group<VertexBase>* > > m_aaGroupInfoVRT;
	Grid::AttachmentAccessor<EdgeBase, Attachment<Group<EdgeBase>* > > m_aaGroupInfoEDG;
	Grid::AttachmentAccessor<Face, Attachment<Group<Face>* > > m_aaGroupInfoFCE;
	Grid::AttachmentAccessor<Volume, Attachment<Group<Volume>* > > m_aaGroupInfoVOL;

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

/*
 * periodic_boundary_manager.h
 *
 *  Created on: 26.11.2012
 *      Author: marscher
 */

#ifndef PERIODIC_IDENTIFIER_H_
#define PERIODIC_IDENTIFIER_H_

#include "lib_grid/grid/grid.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/grid/geometric_base_objects.h"

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
	virtual bool match(Volume*, Volume*) {UG_THROW("not impled, because volume identification is not supported.")}
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
	virtual bool match(VertexBase* v1, VertexBase* v2) {return match_impl(v1, v2);}
	virtual bool match(EdgeBase* e1, EdgeBase* e2) {return match_impl(e1, e2);}
	virtual bool match(Face* f1, Face* f2) {return match_impl(f1, f2);}

	virtual ~ParallelShiftIdentifier() {}
	typedef typename TPosAA::ValueType AttachmentType;
	ParallelShiftIdentifier(TPosAA& aa) : m_aaPos(aa) {}
	void set_shift(AttachmentType& shift) {m_shift = shift; VecScale(m_shift_opposite, m_shift, -1);}
protected:
	AttachmentType m_shift;
	AttachmentType m_shift_opposite;
	TPosAA& m_aaPos;
	template<class TElem> bool match_impl(TElem*, TElem*) const;
};


//template<class TPosAA, int dim> class TransformationBasedIdentifier : public IIdentifier {
//public:
//	TransformationBasedIdentifier(TPosAA& aa) : m_aaPos(aa) {}
//	void setTransformation(MathMatrix<dim,dim>& T) {this->T = T;}
//
//	virtual bool match(VertexBase*, VertexBase*);
//	virtual bool match(EdgeBase*, EdgeBase*);
//	virtual bool match(Face*, Face*);
//protected:
//	MathMatrix<dim,dim> T;
//	TPosAA& m_aaPos;
//};

///
/**
 *
 */
class PeriodicBoundaryManager : public GridObserver {
public:

	/**
	 * A Group instance holds a master of type TElem and several children.
	 */
	template <class TElem, class Container = std::vector<TElem*> > class Group {
	public:
		typedef Container SlaveContainer;
		typedef typename Container::iterator SlaveIterator;

		Group(TElem* m = NULL) : m_master(m) {}

		void add_slave(TElem* e) {
			UG_ASSERT(e, "add_slave: slave not valid");
			UG_ASSERT(e != m_master, "element already master!");
			m_slaves.push_back(e); }

		Container& get_slaves() { return m_slaves; }
		TElem* m_master;

	protected:
		Container m_slaves;
	};

	enum PeriodicStatus {
		P_SLAVE, P_SLAVE_MASTER_UNKNOWN, P_MASTER, P_NOT_PERIODIC
	};

	PeriodicBoundaryManager() : m_pGrid(NULL), m_pIdentifier(NULL) {}
	~PeriodicBoundaryManager();

	// sets grid and inits group info attachment accessors
	void set_grid(Grid* g);

	///
	/**
	 * ATTENTION: method assumes, that e1 and e2 are matching!
	 * identifies element e1 with e2. If element type has sides, it will recursively
	 * ascend to the lowest dimension type (vertex) and identifies them.
	 * @param e1
	 * @param e2
	 * @param i
	 */
	template <class TElem> void identify(TElem* e1, TElem* e2, IIdentifier* i =
			NULL);
	/**
	 * performs matching before identification using given IIdentifier or the local
	 * last used Identifier. Throws exception if no identifier is usable.
	 */
	template <class TElem> void match_and_identify(TElem* e1, TElem* e2,
			IIdentifier* i = NULL);
	template <class TElem> bool is_periodic(TElem* e) const;
	template <class TElem> bool is_slave(TElem*) const;
	template <class TElem> bool is_master(TElem*) const;
	template <class TElem> TElem* master(TElem* e) const;
	template <class TElem> typename Group<TElem>::SlaveContainer* slaves(
			TElem* e) const;
	template <class TElem> void print_identification() const;


	/// grid observation methods
	virtual void grid_to_be_destroyed(Grid* grid);
	virtual void vertex_created(Grid* grid, VertexBase* vrt,
										GeometricObject* pParent = NULL,
										bool replacesParent = false);

	virtual void edge_created(Grid* grid, EdgeBase* e,
								GeometricObject* pParent = NULL,
								bool replacesParent = false);

	virtual void face_created(Grid* grid, Face* f,
								GeometricObject* pParent = NULL,
								bool replacesParent = false);

	virtual void vertex_to_be_erased(Grid* grid, VertexBase* vrt,
									 VertexBase* replacedBy = NULL);

	virtual void edge_to_be_erased(Grid* grid, EdgeBase* e,
									 EdgeBase* replacedBy = NULL);

	virtual void face_to_be_erased(Grid* grid, Face* f,
									 Face* replacedBy = NULL);

	void set_identifier(IIdentifier*);

protected:
	// no copy construction allowed
	PeriodicBoundaryManager(const PeriodicBoundaryManager&)	{}

	/// grid instance we operate on
	MultiGrid* m_pGrid;

	/// identifier used
	std::auto_ptr<IIdentifier> m_pIdentifier;

	/// attachment accessors for Groups
	Grid::AttachmentAccessor<VertexBase, Attachment<Group<VertexBase>* > > m_aaGroupVRT;
	Grid::AttachmentAccessor<EdgeBase, Attachment<Group<EdgeBase>* > > m_aaGroupEDG;
	Grid::AttachmentAccessor<Face, Attachment<Group<Face>* > > m_aaGroupFCE;

	/// attachment accessors for PeriodicStatus
	Grid::AttachmentAccessor<VertexBase, Attachment<PeriodicStatus> > m_aaPeriodicStatusVRT;
	Grid::AttachmentAccessor<EdgeBase, Attachment<PeriodicStatus> > m_aaPeriodicStatusEDG;
	Grid::AttachmentAccessor<Face, Attachment<PeriodicStatus> > m_aaPeriodicStatusFCE;

	/// make element e slave of group g
	template <class TElem> void make_slave(Group<TElem>* g, TElem* e);
	template <class TElem> bool remove_slave(TElem* slave);
	template <class TElem> void merge_groups(Group<TElem>* g0, Group<TElem>* g1);

	// get typed attachment accessor for group attachment
	template <class TElem>
	const Grid::AttachmentAccessor<TElem, Attachment<Group<TElem>* > >&
	get_group_accessor() const;

	template <class TElem>
	Grid::AttachmentAccessor<TElem, Attachment<Group<TElem>* > >&
	get_group_accessor();

	// get typed attachment accessor for periodic status attachment
	template <class TElem>
	const Grid::AttachmentAccessor<TElem, Attachment<PeriodicStatus> >&
	get_periodic_status_accessor() const;

	template <class TElem>
	Grid::AttachmentAccessor<TElem, Attachment<PeriodicStatus> >&
	get_periodic_status_accessor();

	// gets group of given element e
	template <class TElem> Group<TElem>* group(TElem* e) const;
	// set group attachment to element e
	template <class TElem> void set_group(Group<TElem>* g, TElem* e);
	template <class TElem> void remove_group(Group<TElem>* g);

	/// handles creation of element type
	template <class TElem, class TParent>
	void handle_creation(TElem* e, TParent* pParent, bool replacesParent);

	/// handles deletion of element type
	template <class TElem>
	void handle_deletion(TElem* e, TElem* replacedBy);

	/// casts parent pointer to exact type before calling handle_creation
	template <class TElem>
	void handle_creation_cast_wrapper(TElem* e, GeometricObject* parent, bool replacesParent);
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
#include "./periodic_boundary_manager_impl.hpp"

#endif /* PERIODIC_IDENTIFIER_H_ */

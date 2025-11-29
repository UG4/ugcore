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

#ifndef PERIODIC_IDENTIFIER_H_
#define PERIODIC_IDENTIFIER_H_

#include "lib_grid/grid/grid.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/grid/grid_base_objects.h"

#include <set>

namespace ug {

/// Interface to match periodic geometric elements
/**
 *
 */
class IIdentifier {
public:
	virtual ~IIdentifier() = default;
	virtual bool match(Vertex*, Vertex*) = 0;
	virtual bool match(Edge*, Edge*) = 0;
	virtual bool match(Face*, Face*) = 0;
	virtual bool match(Volume*, Volume*) {UG_THROW("not impled, because volume identification is not supported.")}
};

/// This class matches geometric elements which are parallel translated.
/**
 * Usage: class needs to be instantiated with the position attachment used on the Domain.
 * Before using any match methods, the translation vector needs to be set with set_shift()
 * \tparam TPosAA {class needs to be instantiated with the
 * position attachment used on the Domain.}
 */
template <typename TPosAA> class ParallelShiftIdentifier: public IIdentifier {
public:
	bool match(Vertex* v1, Vertex* v2) override {return match_impl(v1, v2);}
	bool match(Edge* e1, Edge* e2) override {return match_impl(e1, e2);}
	bool match(Face* f1, Face* f2) override {return match_impl(f1, f2);}

	~ParallelShiftIdentifier() override = default;

	using AttachmentType = typename TPosAA::ValueType;
	ParallelShiftIdentifier(TPosAA& aa) : m_aaPos(aa) {}
	void set_shift(AttachmentType& shift) {m_shift = shift; VecScale(m_shift_opposite, m_shift, -1);}
protected:
	AttachmentType m_shift;
	AttachmentType m_shift_opposite;
	TPosAA& m_aaPos;
	template<typename TElem> bool match_impl(TElem*, TElem*) const;
};


//template<typename TPosAA, int dim> class TransformationBasedIdentifier : public IIdentifier {
//public:
//	TransformationBasedIdentifier(TPosAA& aa) : m_aaPos(aa) {}
//	void setTransformation(MathMatrix<dim,dim>& T) {this->T = T;}
//
//	virtual bool match(Vertex*, Vertex*);
//	virtual bool match(Edge*, Edge*);
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
	template <typename TElem, typename Container = std::vector<TElem*> >
	class Group
	{
		public:
			using SlaveContainer = Container;
			using SlaveIterator = typename Container::iterator;
			// the set type definitions are used to check for periodicity after identification
			using master_slave_pair = std::pair<TElem*, TElem*>;
			using unique_pairs = std::set<master_slave_pair>;

			Group(TElem* m = nullptr) : m_master(m) {}

			void add_slave(TElem* e) {
				UG_ASSERT(e, "add_slave: slave not valid");
				UG_ASSERT(e != m_master, "element already master!");
				m_slaves.push_back(e); }

			Container& get_slaves() { return m_slaves; }
			TElem* m_master;

		protected:
			Container m_slaves;
	};

	enum PeriodicStatus : byte_t {
		P_SLAVE, P_SLAVE_MASTER_UNKNOWN, P_MASTER, P_NOT_PERIODIC
	};


	PeriodicBoundaryManager();
	~PeriodicBoundaryManager() override;

	// sets grid and inits group info attachment accessors
	void set_grid(Grid* g);
	[[nodiscard]] Grid* get_grid() const;

	// sets the subset handler to use for element lookup
//	void set_subset_handler(ISubsetHandler* sh);

	///
	/**
	 * ATTENTION: method assumes, that e1 and e2 are matching!
	 * identifies element e1 with e2. If element type has sides, it will recursively
	 * ascend to the lowest dimension type (vertex) and identifies them.
	 * @param e1
	 * @param e2
	 * @param i identifier to use to perform geometrical matching of elements
	 */
	template <typename TElem> void identify(TElem* e1, TElem* e2, IIdentifier& i);

	template <typename TElem> bool is_periodic(TElem* e) const;
	template <typename TElem> bool is_slave(TElem*) const;
	template <typename TElem> bool is_master(TElem*) const;
	template <typename TElem> TElem* master(TElem* e) const;
	template <typename TElem> typename Group<TElem>::SlaveContainer* slaves( TElem* e) const;
	template <typename TElem> void print_identification() const;


	/// grid observation methods
	void grid_to_be_destroyed(Grid* grid) override;

	void vertex_created(Grid* grid, Vertex* vrt,
	                    GridObject* pParent = nullptr,
	                    bool replacesParent = false) override;

	void edge_created(Grid* grid, Edge* e,
	                  GridObject* pParent = nullptr,
	                  bool replacesParent = false) override;

	void face_created(Grid* grid, Face* f,
	                  GridObject* pParent = nullptr,
	                  bool replacesParent = false) override;

	void vertex_to_be_erased(Grid* grid, Vertex* vrt,
	                         Vertex* replacedBy = nullptr) override;

	void edge_to_be_erased(Grid* grid, Edge* e,
	                       Edge* replacedBy = nullptr) override;

	void face_to_be_erased(Grid* grid, Face* f,
	                       Face* replacedBy = nullptr) override;

	/// checks that all elements of given gocs are periodic (called after identification)
	bool check_periodicity(const GridObjectCollection& goc1,
							  const GridObjectCollection& goc2,
							  ISubsetHandler* sh = nullptr);

	/**	makes sure that all master/slave identifications are correct and that
	 * no unconnected masters or slaves exist.
	 * This method iterates over all elements of the grid and is thus rather expensive.*/
	void validity_check();

	/**
	 * sets the identifier instance to use for subset index si
	 * @param i identifier pointer
	 * @param si
	 */
//	void set_identifier(SmartPtr<IIdentifier> i, size_t si);

protected:
	// no copy construction allowed
	PeriodicBoundaryManager(const PeriodicBoundaryManager&)	= delete;

	/// grid instance we operate on
	MultiGrid* m_pGrid;

	/// store subset handler of domain to lookup element subset ids
//	ISubsetHandler* m_pSH;

	/// identifier mapping to subset indices
	/**
	 * stores IIdentifier pointers to subset index position to look them up
	 * in constant time.
	 */
//	std::vector<SmartPtr<IIdentifier> > m_vIdentifier;

	/// attachment accessors for Groups
	Grid::AttachmentAccessor<Vertex, Attachment<Group<Vertex>* > > m_aaGroupVRT;
	Grid::AttachmentAccessor<Edge, Attachment<Group<Edge>* > > m_aaGroupEDG;
	Grid::AttachmentAccessor<Face, Attachment<Group<Face>* > > m_aaGroupFCE;

	/// attachment accessors for PeriodicStatus
	Grid::AttachmentAccessor<Vertex, Attachment<PeriodicStatus> > m_aaPeriodicStatusVRT;
	Grid::AttachmentAccessor<Edge, Attachment<PeriodicStatus> > m_aaPeriodicStatusEDG;
	Grid::AttachmentAccessor<Face, Attachment<PeriodicStatus> > m_aaPeriodicStatusFCE;

	/// make element e slave of group g
	template <typename TElem> void make_slave(Group<TElem>* g, TElem* e);
	template <typename TElem> bool remove_slave(TElem* slave);

	// make e a master of g
	template <typename TElem> void make_master(Group<TElem>* g, TElem* e);

	template <typename TElem> void merge_groups(Group<TElem>* g0, Group<TElem>* g1);

	// get typed attachment accessor for group attachment
	template <typename TElem>
	const Grid::AttachmentAccessor<TElem, Attachment<Group<TElem>* > >&
	get_group_accessor() const;

	template <typename TElem>
	Grid::AttachmentAccessor<TElem, Attachment<Group<TElem>* > >&
	get_group_accessor();

	// get typed attachment accessor for periodic status attachment
	template <typename TElem>
	const Grid::AttachmentAccessor<TElem, Attachment<PeriodicStatus> >&
	get_periodic_status_accessor() const;

	template <typename TElem>
	Grid::AttachmentAccessor<TElem, Attachment<PeriodicStatus> >&
	get_periodic_status_accessor();

	// gets group of given element e
	template <typename TElem> Group<TElem>* group(TElem* e) const;
	// set group attachment to element e
	template <typename TElem> void set_group(Group<TElem>* g, TElem* e);
	template <typename TElem> void remove_group(Group<TElem>* g);

	///	replaces all group occurrences of pParent by the specified elem
	template <typename TElem>
	void replace_parent(TElem* e, TElem* pParent);

	/// handles creation of element type
	template <typename TElem, typename TParent>
	void handle_creation(TElem* e, TParent* pParent);

	/// handles deletion of element type
	template <typename TElem>
	void handle_deletion(TElem* e, TElem* replacedBy);

	/// casts parent pointer to exact type before calling handle_creation
	template <typename TElem>
	void handle_creation_cast_wrapper(TElem* e, GridObject* parent, bool replacesParent);

	template <typename TElem, typename TIterator>
	void check_elements_periodicity(
			TIterator begin,
			TIterator end,
			typename Group<TElem>::unique_pairs& s,
			ISubsetHandler* sh);

	template <typename TElem>
	void validity_check();

};

///	Accesses attachments with consideration to periodic boundaries
/**	Slave element redirects to related master element
 *  also works in case of no periodic boundary conditions
 *  (but is then probably slower than standard attachment accessor)
 */
template <typename TElem, typename TAttachment>
class PeriodicAttachmentAccessor
{
	public:
		using ValueType = typename TAttachment::ValueType;
		using RefType = typename attachment_value_traits<ValueType>::reference;
		using ConstRefType = typename attachment_value_traits<ValueType>::const_reference;

		PeriodicAttachmentAccessor() : m_pbm(nullptr)	{}

		PeriodicAttachmentAccessor(Grid& g, TAttachment& a) : m_pbm(nullptr)
		{
			access(g, a);
		}

		explicit PeriodicAttachmentAccessor(PeriodicBoundaryManager& pbm) : m_pbm(&pbm) {}

		bool access(Grid& g, TAttachment& a)
		{
			if (!(m_pbm)) m_pbm = g.periodic_boundary_manager();
			return m_aa.access(g, a);
		}

		RefType operator [] (TElem* e)	{
			if(m_pbm && m_pbm->is_slave(e))
				return m_aa[m_pbm->master(e)];
			return m_aa[e];
		}

		ConstRefType operator [] (TElem* e) const	{
			if(m_pbm && m_pbm->is_slave(e))
				return m_aa[m_pbm->master(e)];
			return m_aa[e];
		}

	private:
		Grid::AttachmentAccessor<TElem, TAttachment>	m_aa;
		PeriodicBoundaryManager* m_pbm;
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
template <typename TDomain>
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
template <typename TDomain>
void IdentifySubsets(TDomain& dom, const char* sName1, const char* sName2);

} // end of namespace ug

// include implementation
#include "./periodic_boundary_manager_impl.hpp"

#endif
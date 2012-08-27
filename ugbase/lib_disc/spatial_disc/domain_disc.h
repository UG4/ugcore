/*
 * domain_disc.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__

// other ug4 modules
#include "common/common.h"
#include "common/util/string_util.h"

// library intern headers
#include "subset_assemble_util.h"
#include "domain_disc_interface.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_assemble_util.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "disc_item.h"
#include "lib_disc/spatial_disc/domain_disc_base.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// domain discretization implementing the interface
/**
 * This class is an implementation of the IDomainDiscretization interface. It
 * is designed to simply group several discreizations on different subsets and
 * perform element based assemblings and costraints in the same order.
 */
template <typename TDomain, typename TAlgebra>
class DomainDiscretization
 :	public DomainDiscBase<TDomain,TAlgebra,
  							DomainDiscretization<TDomain, TAlgebra> >
{
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of approximation space
		typedef ApproximationSpace<TDomain>	approx_space_type;

	public:
	///	default Constructor
		DomainDiscretization(SmartPtr<approx_space_type> pApproxSpace) :
			m_spApproxSpace(pApproxSpace), m_bForceRegGrid(false),
			m_ConstraintTypesEnabled(CT_ALL), m_ElemTypesEnabled(EDT_ALL),
			m_pBoolMarker(NULL)
		{
			this->set_approximation_space(pApproxSpace);
		};

	///////////////////////////
	// Time independent part
	///////////////////////////

	/// \copydoc IAssemble::assemble_jacobian()
		template <typename TDD>
		void assemble_jacobian(matrix_type& J, const vector_type& u, ConstSmartPtr<TDD> dd);

	/// \copydoc IAssemble::assemble_defect()
		template <typename TDD>
		void assemble_defect(vector_type& d, const vector_type& u, ConstSmartPtr<TDD> dd);

	/// \copydoc IAssemble::assemble_linear()
		template <typename TDD>
		void assemble_linear(matrix_type& A, vector_type& b, ConstSmartPtr<TDD> dd);

	/// assembles the stiffness matrix
		template <typename TDD>
		void assemble_rhs(vector_type& rhs, const vector_type& u, ConstSmartPtr<TDD> dd);

	/// \copydoc IAssemble::assemble_rhs()
		template <typename TDD>
		void assemble_rhs(vector_type& b, ConstSmartPtr<TDD> dd);

	/// \copydoc IAssemble::adjust_solution()
		template <typename TDD>
		void adjust_solution(vector_type& u, ConstSmartPtr<TDD> dd);

	///////////////////////
	// Time dependent part
	///////////////////////

	/// \copydoc IDomainDiscretization::prepare_timestep()
		template <typename TDD>
		void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                      ConstSmartPtr<TDD> dd);

	/// \copydoc IDomainDiscretization::assemble_jacobian()
		template <typename TDD>
		void assemble_jacobian(matrix_type& J,
		                       ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                       const number s_a0,
		                       ConstSmartPtr<TDD> dd);

	/// \copydoc IDomainDiscretization::assemble_defect()
		template <typename TDD>
		void assemble_defect(vector_type& d,
		                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff,
		                     ConstSmartPtr<TDD> dd);

	/// \copydoc IDomainDiscretization::assemble_linear()
		template <typename TDD>
		void assemble_linear(matrix_type& A, vector_type& b,
		                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff,
		                     ConstSmartPtr<TDD> dd);

	/// \copydoc IDomainDiscretization::assemble_rhs()
		template <typename TDD>
		void assemble_rhs(	 vector_type& b,
							 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
							 const std::vector<number>& vScaleMass,
							 const std::vector<number>& vScaleStiff,
							 ConstSmartPtr<TDD> dd);

	/// \copydoc IDomainDiscretization::adjust_solution()
		template <typename TDD>
		void adjust_solution(vector_type& u, number time,
		                       ConstSmartPtr<TDD> dd);

	/// \copydoc IDomainDiscretization::finish_timestep()
		template <typename TDD>
		void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
							 ConstSmartPtr<TDD> dd);

	///////////////////////////
	// Mass and Stiffness Matrix
	///////////////////////////

	/// assembles the mass matrix
		template <typename TDD>
		void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                          ConstSmartPtr<TDD> dd);

	/// assembles the stiffness matrix
		template <typename TDD>
		void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                               ConstSmartPtr<TDD> dd);

	public:
	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce) {m_bForceRegGrid = bForce;}

	///	returns if constraints enabled
		int constraints_enabled() const {return m_ConstraintTypesEnabled;}

	///	enables constraints
		void enable_constraints(int bEnableTypes) {m_ConstraintTypesEnabled = bEnableTypes;}

	///	returns type of boundary elem discs enabled
		virtual int elem_discs_enabled() const {return m_ElemTypesEnabled;}

	///	enables boundary elem discs
		virtual void enable_elem_discs(int bEnableTypes) {m_ElemTypesEnabled = bEnableTypes;}

	///	sets a selector to exlude elements from assembling
	/**
	 * This methods sets a selector. Only elements that are selected will be
	 * assembled during assembling process. If no selector is set, this
	 * corresponds to a selector where all elements have been selected.
	 *
	 * \param[in]	sel		Selector
	 */
	virtual void set_selector(BoolMarker* sel = NULL){m_pBoolMarker = sel;}

	public:
	/// adds an element discretization to the assembling process
	/**
	 * This function adds an Element Discretization to the assembling. During
	 * the assembling process the elem disc will assemble a given group of
	 * subsets and add the output to the passed vectors and matrices.
	 * For each Elem Disc one loop over all elements of the subset will be
	 * performed.
	 *
	 * \param[in] 	elem		Element Discretization to be added
	 */
		void add(SmartPtr<IDomainElemDisc<TDomain> > elem)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vDomainElemDisc.size(); ++i)
				if(m_vDomainElemDisc[i] == elem)
					return;

		//	set approximation space
			elem->set_approximation_space(m_spApproxSpace);

		//	add it
			m_vDomainElemDisc.push_back(elem);
		}

	/// adds a constraint to the assembling process
	/**
	 * This function adds a IConstraint to the assembling. The constraint is
	 * called when all element-wise assembling have been performed.
	 *
	 * \param[in] 	pp		Constraint to be added
	 */
		void add(SmartPtr<IDomainConstraint<TDomain, TAlgebra> > pp)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i] == pp)
					return;

		//	add constraint
			m_vConstraint.push_back(pp);
		}

	/// adds a disc item to the assembling process
	/**
	 * This function adds a IDiscretizationItem to the assembling. The contained
	 * elem discs and constraints are sorted into the lists
	 *
	 * \param[in] 	di		Disc Item
	 */
		void add(SmartPtr<IDiscretizationItem<TDomain, TAlgebra> > di)
		{
		//	add elem discs
			for(size_t i = 0; i < di->num_elem_disc(); ++i)
				add(di->elem_disc(i));

		//	add constraints
			for(size_t i = 0; i < di->num_constraint(); ++i)
				add(di->constraint(i));
		}

	protected:
	///	returns number of registered dirichlet constraints
		virtual size_t num_dirichlet_constraints() const
		{
			size_t cnt = 0;
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i]->type() & CT_DIRICHLET)
					cnt++;

			return cnt;
		}

	///	returns the i'th dirichlet constraint
		virtual SmartPtr<IConstraint<TAlgebra> > dirichlet_constraint(size_t nr)
		{
			size_t cnt = 0;
			for(size_t i = 0; i < m_vConstraint.size(); ++i)
				if(m_vConstraint[i]->type() & CT_DIRICHLET)
				{
					if(cnt == nr) return m_vConstraint[i];
					cnt++;
				}
			UG_THROW("DomainDisc: Dirichlet Constraint "<<nr<<" not found.");
		}

	///	set the approximation space in the elem discs and extract IElemDiscs
		void update_elem_discs();
		void update_constraints();
		void update_disc_items();

	protected:
	///	vector holding all registered elem discs
		std::vector<SmartPtr<IDomainElemDisc<domain_type> > > m_vDomainElemDisc;

	///	vector holding all registered elem discs
		std::vector<IElemDisc*> m_vElemDisc;

	//	vector holding all registered constraints
		std::vector<SmartPtr<IDomainConstraint<TDomain, TAlgebra> > > m_vConstraint;

	///	current approximation space
		SmartPtr<approx_space_type> m_spApproxSpace;

	/// forces the assembling to regard the grid as regular
		bool m_bForceRegGrid;

	///	enables the constraints
		int m_ConstraintTypesEnabled;

	///	enables the constraints
		int m_ElemTypesEnabled;

	///	selector used to skip elements
		BoolMarker* m_pBoolMarker;
};

/// @}

} // end namespace ug

// include documentation
#include "domain_disc_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__ */

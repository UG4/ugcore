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
#include "lib_disc/function_spaces/grid_function.h"

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
			m_pBoolMarker(NULL), m_pSelector(NULL)
		{
			this->set_approximation_space(pApproxSpace);
		};

		virtual ~DomainDiscretization() {};

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

	///	wrapper for GridFunction
	/// \{
		template <typename TDD>
		void assemble_jacobian(matrix_type& J, GridFunction<TDomain, TDD, TAlgebra>& u)
			{assemble_jacobian<TDD>(J, u, u.dof_distribution());}

		template <typename TDD>
		void assemble_defect(vector_type& d, GridFunction<TDomain, TDD, TAlgebra>& u)
			{assemble_defect<TDD>(d, u, u.dof_distribution());}

		template <typename TDD>
		void assemble_linear(matrix_type& A, GridFunction<TDomain, TDD, TAlgebra>& rhs)
			{assemble_linear<TDD>(A, rhs, rhs.dof_distribution());}

		template <typename TDD>
		void assemble_rhs(vector_type& rhs, GridFunction<TDomain, TDD, TAlgebra>& u)
			{assemble_rhs<TDD>(rhs, u, u.dof_distribution());}

		template <typename TDD>
		void assemble_rhs(GridFunction<TDomain, TDD, TAlgebra>& b)
			{assemble_rhs<TDD>(b, b.dof_distribution());}

		template <typename TDD>
		void adjust_solution(GridFunction<TDomain, TDD, TAlgebra>& u)
			{adjust_solution<TDD>(u, u.dof_distribution());}
	/// \}

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

	///	sets a marker to exlude elements from assembling
	/**
	 * This methods sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
	virtual void set_marker(BoolMarker* mark = NULL){m_pBoolMarker = mark;}

	///	sets a selector of elements for assembling
	/**
	 * This methods sets an element list. Only elements of this list will be
	 * assembled during assembling process. Especially the list defines the begin
	 * and end of the element-iterator in the element assembling-loop.
	 * If no element list is set, this corresponds to a assembling where the loop is
	 * carried out over all elements of a subset.
	 *
	 * \param[in]	sel		Selector
	 */
	virtual void set_selector(Selector* sel = NULL){m_pSelector = sel;}

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

	///	returns number of registered constraints
		virtual size_t num_constraints() const {return m_vConstraint.size();}

	///	returns the i'th constraint
		virtual SmartPtr<IConstraint<TAlgebra> > constraint(size_t i) {return m_vConstraint[i];}

	protected:
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

	///	marker used to skip elements
		BoolMarker* m_pBoolMarker;
	///	selector used to set a list of elements for the assembling
		Selector* 	m_pSelector;
};

/// @}

} // end namespace ug

// include documentation
#include "domain_disc_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__ */

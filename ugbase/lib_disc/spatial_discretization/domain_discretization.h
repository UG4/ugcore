/*
 * domain_discretization.h
 *
 *  Created on: 29.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__
#define __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__

// other ug4 modules
#include "common/common.h"
#include "common/util/string_util.h"

// library intern headers
#include "subset_assemble_util.h"
#include "domain_discretization_interface.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/spatial_discretization/elem_disc/elem_disc_assemble_util.h"
#include "lib_disc/spatial_discretization/constraints/constraint_interface.h"
#include "disc_item.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// domain discretization implementing the interface
/**
 * This class is an implementation of the IDomainDiscretization interface. It
 * is designed to simply group several discreizations on different subsets and
 * perform element based assemblings and costraints in the same order.
 */
template <	typename TDomain,
			typename TDoFDistribution,
			typename TAlgebra>
class DomainDiscretization :
	public IDomainDiscretization<TDoFDistribution, TAlgebra>
{
	public:
	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of DoF Distribution
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	///	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	///	Type of approximation space
		typedef ApproximationSpace<TDomain, TDoFDistribution, TAlgebra>
					approx_space_type;

	public:
	///	default Constructor
		DomainDiscretization() :
			m_pApproxSpace(NULL), m_bForceRegGrid(false), m_pSelector(NULL) {};

	///	sets the approximation space to use
		void set_approximation_space(approx_space_type& pApproxSpace)
		{
			m_pApproxSpace = & pApproxSpace;
		}

	///////////////////////////
	// Time independent part
	///////////////////////////

	/// \copydoc IAssemble::assemble_jacobian()
		bool assemble_jacobian(matrix_type& J, const vector_type& u,
	                           const dof_distribution_type& dd);

		bool assemble_jacobian(matrix_type& J, const vector_type& u)
			{return assemble_jacobian(J, u, get_surface_dd());}

	/// \copydoc IAssemble::assemble_defect()
		bool assemble_defect(vector_type& d, const vector_type& u,
		                     const dof_distribution_type& dd);

		bool assemble_defect(vector_type& d, const vector_type& u)
			{return assemble_defect(d, u, get_surface_dd());}

	/// \copydoc IAssemble::assemble_linear()
		bool assemble_linear(matrix_type& A, vector_type& b,
		                     const vector_type& u,
		                     const dof_distribution_type& dd);

		bool assemble_linear(matrix_type& A, vector_type& b, const vector_type& u)
			{return assemble_linear(A, b, u, get_surface_dd());}

	/// \copydoc IAssemble::assemble_solution()
		bool assemble_solution(vector_type& u,
		                       const dof_distribution_type& dd);

		bool assemble_solution(vector_type& u)
			{return assemble_solution(u, get_surface_dd());}

	///////////////////////
	// Time dependent part
	///////////////////////

	/// \copydoc IDomainDiscretization::assemble_jacobian()
		bool assemble_jacobian(matrix_type& J,
		                       const VectorTimeSeries<vector_type>& vSol,
		                       const number s_a0,
		                       const dof_distribution_type& dd);

		bool assemble_jacobian(matrix_type& J,
		                       const VectorTimeSeries<vector_type>& vSol,
		                       const number s_a0)
		{return assemble_jacobian(J, vSol, s_a0, get_surface_dd());}

	/// \copydoc IDomainDiscretization::assemble_defect()
		bool assemble_defect(vector_type& d,
		                     const VectorTimeSeries<vector_type>& vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff,
		                     const dof_distribution_type& dd);

		bool assemble_defect(vector_type& d,
		                     const VectorTimeSeries<vector_type>& vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff)
		{return assemble_defect(d, vSol, vScaleMass, vScaleStiff, get_surface_dd());}

	/// \copydoc IDomainDiscretization::assemble_linear()
		bool assemble_linear(matrix_type& A, vector_type& b,
		                     const VectorTimeSeries<vector_type>& vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff,
		                     const dof_distribution_type& dd);

		bool assemble_linear(matrix_type& A, vector_type& b,
		                     const VectorTimeSeries<vector_type>& vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff)
		{return assemble_linear(A, b, vSol, vScaleMass, vScaleStiff, get_surface_dd());}

	/// \copydoc IDomainDiscretization::assemble_solution()
		bool assemble_solution(vector_type& u, number time,
		                       const dof_distribution_type& dd);

		bool assemble_solution(vector_type& u, number time)
		{return assemble_solution(u, time, get_surface_dd());}

	///////////////////////////
	// Mass and Stiffness Matrix
	///////////////////////////

	/// assembles the mass matrix
		bool assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                          const dof_distribution_type& dd);

		bool assemble_mass_matrix(matrix_type& M, const vector_type& u)
		{return assemble_mass_matrix(M, u, get_surface_dd());}

	/// assembles the stiffness matrix
		bool assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                               const dof_distribution_type& dd);

		bool assemble_stiffness_matrix(matrix_type& A, const vector_type& u)
		{return assemble_stiffness_matrix(A, u, get_surface_dd());}

	/// assembles the stiffness matrix
		bool assemble_rhs(vector_type& rhs, const vector_type& u,
		                  const dof_distribution_type& dd);

		bool assemble_rhs(vector_type& rhs, const vector_type& u)
		{return assemble_rhs(rhs, u, get_surface_dd());}

	public:
	/// forces the assembling to consider the grid as regular
		virtual void force_regular_grid(bool bForce) {m_bForceRegGrid = bForce;}

	///	sets a selector to exlude elements from assembling
	/**
	 * This methods sets a selector. Only elements that are selected will be
	 * assembled during assembling process. If no selector is set, this
	 * corresponds to a selector where all elements have been selected.
	 *
	 * \param[in]	sel		Selector
	 */
	virtual void set_selector(ISelector* sel = NULL){m_pSelector = sel;}

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
		bool add(IDomainElemDisc<domain_type>& elem)
		{
		//	check that not already registered
			for(size_t i = 0; i < m_vDomainElemDisc.size(); ++i)
				if(m_vDomainElemDisc[i] == &elem)
					return true;

		//	add it
			m_vDomainElemDisc.push_back(&elem);
			return true;
		}

	/// adds a constraint to the assembling process
	/**
	 * This function adds a IConstraint to the assembling. The constraint is
	 * called when all element-wise assembling have been performed.
	 *
	 * \param[in] 	pp		Constraint to be added
	 */
		bool add(IConstraint<TDoFDistribution, TAlgebra>& pp)
		{
		// 	get type of constraint
			const int type = pp.type();

		//	check that not already registered
			for(size_t i = 0; i < m_vvConstraints[type].size(); ++i)
				if(m_vvConstraints[type][i] == &pp)
					return true;

		//	add constraint
			m_vvConstraints[type].push_back(&pp);
			return true;
		}

	/// adds a disc item to the assembling process
	/**
	 * This function adds a IDiscretizationItem to the assembling. The contained
	 * elem discs and constraints are sorted into the lists
	 *
	 * \param[in] 	di		Disc Item
	 */
		bool add(IDiscretizationItem<TDomain, TDoFDistribution, TAlgebra>& di)
		{
			bool bRet = true;

		//	add elem discs
			for(size_t i = 0; i < di.num_elem_disc(); ++i)
				bRet &= add(*di.get_elem_disc(i));

		//	add constraints
			for(size_t i = 0; i < di.num_constraint(); ++i)
				bRet &= add(*di.get_constraint(i));

			return bRet;
		}

	protected:
	///	returns number of registered dirichlet constraints
		virtual size_t num_dirichlet_constraints() const
		{
			return m_vvConstraints[CT_DIRICHLET].size();
		}

	///	returns the i'th dirichlet constraint
		virtual IConstraint<TDoFDistribution, TAlgebra>* get_dirichlet_constraint(size_t i)
		{
			return m_vvConstraints[CT_DIRICHLET].at(i);
		}

	///	set the approximation space in the elem discs and extract IElemDiscs
		bool update_elem_discs();

	///	returns the surface dof distribution
		dof_distribution_type& get_surface_dd();

	protected:
	///	vector holding all registered elem discs
		std::vector<IDomainElemDisc<domain_type>*> m_vDomainElemDisc;

	///	vector holding all registered elem discs
		std::vector<IElemDisc*> m_vElemDisc;

	//	vector holding all registered constraints
		std::vector<IConstraint<TDoFDistribution, TAlgebra>*> m_vvConstraints[NUM_CONSTRAINT_TYPES];

	///	current approximation space
		approx_space_type* m_pApproxSpace;

	/// forces the assembling to regard the grid as regular
		bool m_bForceRegGrid;

	///	selector used to skip elements
		ISelector* m_pSelector;
};

/// @}

} // end namespace ug

// inlcude documentation
#include "domain_discretization_impl.h"

#endif /* __H__UG_LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DOMAIN_DISCRETIZATION__ */

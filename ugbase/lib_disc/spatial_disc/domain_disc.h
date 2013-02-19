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
#include "lib_disc/spatial_disc/ass_adapter.h"

namespace ug {

/// \ingroup lib_disc_domain_assemble
/// @{

/// domain discretization implementing the interface
/**
 * This class is an implementation of the IDomainDiscretization interface. It
 * is designed to simply group several discretizations on different subsets and
 * perform element based assemblings and constraints in the same order.
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
			m_AssAdapter()
		{
			this->set_approximation_space(pApproxSpace);
			
			m_AssAdapter.pBoolMarker = NULL;
			m_AssAdapter.pSelector = NULL;
			m_AssAdapter.assIndex.index_set = false;
		};

		virtual ~DomainDiscretization() {};

	///////////////////////////
	// Time independent part
	///////////////////////////

	/// \copydoc IAssemble::assemble_jacobian()
		void assemble_jacobian(matrix_type& J, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IAssemble::assemble_defect()
		void assemble_defect(vector_type& d, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IAssemble::assemble_linear()
		void assemble_linear(matrix_type& A, vector_type& b, ConstSmartPtr<DoFDistribution> dd);

	/// assembles the stiffness matrix
		void assemble_rhs(vector_type& rhs, const vector_type& u, ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IAssemble::assemble_rhs()
		void assemble_rhs(vector_type& b, ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IAssemble::adjust_solution()
		void adjust_solution(vector_type& u, ConstSmartPtr<DoFDistribution> dd);

	///	wrapper for GridFunction
	/// \{
		void assemble_jacobian(matrix_type& J, GridFunction<TDomain, TAlgebra>& u)
			{assemble_jacobian(J, u, u.dof_distribution());}

		void assemble_defect(vector_type& d, GridFunction<TDomain, TAlgebra>& u)
			{assemble_defect(d, u, u.dof_distribution());}

		void assemble_linear(matrix_type& A, GridFunction<TDomain, TAlgebra>& rhs)
			{assemble_linear(A, rhs, rhs.dof_distribution());}

		void assemble_rhs(vector_type& rhs, GridFunction<TDomain, TAlgebra>& u)
			{assemble_rhs(rhs, u, u.dof_distribution());}

		void assemble_rhs(GridFunction<TDomain, TAlgebra>& b)
			{assemble_rhs(b, b.dof_distribution());}

		void adjust_solution(GridFunction<TDomain, TAlgebra>& u)
			{adjust_solution(u, u.dof_distribution());}
	/// \}

	///////////////////////
	// Time dependent part
	///////////////////////

	/// \copydoc IDomainDiscretization::prepare_timestep()
		void prepare_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                      ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IDomainDiscretization::assemble_jacobian()
		void assemble_jacobian(matrix_type& J,
		                       ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                       const number s_a0,
		                       ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IDomainDiscretization::assemble_defect()
		void assemble_defect(vector_type& d,
		                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff,
		                     ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IDomainDiscretization::assemble_linear()
		void assemble_linear(matrix_type& A, vector_type& b,
		                     ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
		                     const std::vector<number>& vScaleMass,
		                     const std::vector<number>& vScaleStiff,
		                     ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IDomainDiscretization::assemble_rhs()
		void assemble_rhs(	 vector_type& b,
							 ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
							 const std::vector<number>& vScaleMass,
							 const std::vector<number>& vScaleStiff,
							 ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IDomainDiscretization::adjust_solution()
		void adjust_solution(vector_type& u, number time,
		                       ConstSmartPtr<DoFDistribution> dd);

	/// \copydoc IDomainDiscretization::finish_timestep()
		void finish_timestep(ConstSmartPtr<VectorTimeSeries<vector_type> > vSol,
							 ConstSmartPtr<DoFDistribution> dd);

	///////////////////////////
	// Mass and Stiffness Matrix
	///////////////////////////

	/// assembles the mass matrix
		void assemble_mass_matrix(matrix_type& M, const vector_type& u,
		                          ConstSmartPtr<DoFDistribution> dd);

	/// assembles the stiffness matrix
		void assemble_stiffness_matrix(matrix_type& A, const vector_type& u,
		                               ConstSmartPtr<DoFDistribution> dd);

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

	///	sets a marker to exclude elements from assembling
	/**
	 * This methods sets a marker. Only elements that are marked will be
	 * assembled during assembling process. If no marker is set, this
	 * corresponds to a marker where all elements have been marked.
	 *
	 * \param[in]	mark	BoolMarker
	 */
	virtual void set_marker(BoolMarker* mark = NULL){m_AssAdapter.pBoolMarker = mark;}

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
	virtual void set_selector(Selector* sel = NULL){
		m_AssAdapter.pSelector = sel;}

	///	sets an index for which the assembling should be carried out
	/**
	 * This methods sets a boolean if an index-wise assemble routine should be used.
	 * This proceeding is e.g. useful for a nonlinear Gauss-Seidel or nonlinear
	 * Jacobi solver. The specific index is passed to the domain discretization.
	 *
	 * \param[in]	ind			size_t
	 * \param[in]	index_set	bool
	 */
	virtual void ass_index(){ ass_index(0.0, false);}
	virtual void ass_index(size_t ind, bool index_set = true)
	{
		m_AssAdapter.assIndex.index = ind; m_AssAdapter.assIndex.index_set = index_set;
	}
	
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
		
	///	this object provides tools to adapt the assemble routine
		AssAdapter m_AssAdapter;
};

/// @}

} // end namespace ug

// include documentation
#include "domain_disc_impl.h"

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DOMAIN_DISC__ */

/*
 * assembled_linear_operator.h
 *
 *  Created on: ..
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__

#include "lib_algebra/operator/operator_interface.h"

#ifdef UG_PARALLEL
#include "lib_discretization/parallelization/parallelization_util.h"
#endif

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
class AssembledLinearOperator :
	public virtual IMatrixOperator<	typename TAlgebra::vector_type,
									typename TAlgebra::vector_type,
									typename TAlgebra::matrix_type>
{
	public:
	// 	Type of algebra
		typedef TAlgebra algebra_type;

	//	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	//	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	default Constructor
		AssembledLinearOperator() :
			m_bInit(false), m_bAssembleRhs(false),
			m_pAss(NULL), m_pDoFDistribution(NULL)
			{};

	///	Constructor
		AssembledLinearOperator(IAssemble<TDoFDistribution, algebra_type>& ass,
		                        bool assemble_rhs = false)
		:	m_bInit(false), m_bAssembleRhs(assemble_rhs),
			m_pAss(&ass), m_pDoFDistribution(NULL)
		{};

	///	sets the discretization to be used
		void set_discretization(IAssemble<TDoFDistribution, algebra_type>& ass)
			{m_pAss = &ass;}

	///	flags, if rhs should be assembled as well
		void export_rhs(bool assemble_rhs) {m_bAssembleRhs = assemble_rhs;}

	///	sets the dof distribution used for assembling
		bool set_dof_distribution(const IDoFDistribution<TDoFDistribution>& dofDistr)
			{m_pDoFDistribution = &dofDistr; return true;}

	///	returns the dof distribution
		const IDoFDistribution<TDoFDistribution>* get_dof_distribution()
				{return m_pDoFDistribution;}

	///	initializes the operator that may depend on the current solution
		virtual bool init(const vector_type& u);

	///	initialize the operator
		virtual bool init();

	///	compute d = J(u)*c (here, J(u) is a Matrix)
		virtual bool apply(vector_type& d, const vector_type& c);

	///	Compute d := d - J(u)*c
		virtual bool apply_sub(vector_type& d, const vector_type& c);

	///	Export matrix
		virtual matrix_type& get_matrix() {return m_J;}

	///	Export assembled rhs
		const vector_type& get_rhs() const {return m_rhs;}

	///	Set Dirichlet values
		bool set_dirichlet_values(vector_type& u);

	/// forces the disc to consider the grid as regular
		void force_regular_grid(bool bForce)
		{
			if(m_pAss != NULL)
				m_pAss->force_regular_grid(bForce);
		}

	///	Destructor
		virtual ~AssembledLinearOperator()
		{
			m_J.destroy();
			m_rhs.destroy();
		};

	protected:
	// 	init flag
		bool m_bInit;

	// 	assemble rhs flag
		bool m_bAssembleRhs;

	// 	assembling procedure
		IAssemble<TDoFDistribution, algebra_type>* m_pAss;

	// 	DoF Distribution used
		const IDoFDistribution<TDoFDistribution>* m_pDoFDistribution;

	// 	matrix storage
		matrix_type m_J;

	// 	vector storage
		vector_type m_rhs;
};

} // namespace ug

// include implementation
#include "assembled_linear_operator_impl.h"

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR__ */

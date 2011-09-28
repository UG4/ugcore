/*
 * assembled_non_linear_operator.h
 *
 *  Created on: ..
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__

#include "lib_algebra/operator/operator_interface.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
class AssembledOperator : public IOperator<	typename TAlgebra::vector_type,
											typename TAlgebra::vector_type>
{
public:
	/// Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Type of DoFDistribution
		typedef TDoFDistribution dof_distribution_type;

	public:
	///	default constructor
		AssembledOperator() :
			m_bInit(false), m_pAss(NULL), m_pDoFDistribution(NULL)
		{};

	///	constructor
		AssembledOperator(IAssemble<dof_distribution_type, algebra_type>& ass) :
			m_bInit(false), m_pAss(&ass), m_pDoFDistribution(NULL)
		{};

	///	sets discretization for assembling
		void set_discretization(IAssemble<TDoFDistribution, algebra_type>& ass) {m_pAss = &ass;}

	///	sets dof distribution
		bool set_dof_distribution(const IDoFDistribution<TDoFDistribution>& dofDistr)
		{
			m_pDoFDistribution = &dofDistr;
			return true;
		}

	///	returns dof distribution
		const IDoFDistribution<TDoFDistribution>* get_dof_distribution()
				{return m_pDoFDistribution;}

	///	Init
		virtual bool init();

	///	Prepare for apply
		virtual bool prepare(vector_type& d, vector_type& u);

	/// Compute d = L(u)
		virtual bool apply(vector_type& d, const vector_type& u);

	/// return assembling
		IAssemble<TDoFDistribution, algebra_type>* get_assemble(){return m_pAss;}

	protected:
		// init flag
		bool m_bInit;

		// assembling procedure
		IAssemble<dof_distribution_type, algebra_type>* m_pAss;

		// DoF Distribution used
		const IDoFDistribution<TDoFDistribution>* m_pDoFDistribution;
};

} // end namepace ug

// include implementation
#include "assembled_non_linear_operator_impl.h"

#endif /*__H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR__*/

/*
 * assembled_linear_operator_impl.h
 *
 *  Created on: 24.02.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__
#define __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__

#include "assembled_linear_operator.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
init(const vector_type& u)
{
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::init: "
				"DoF Distribution not set.\n");
		return false;
	}

	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::init:"
				" Assembling rountine not set.\n");
		return false;
	}

//	get number of dofs
	const size_t numDoFs = m_pDoFDistribution->num_dofs();

//	resize matrix and set to zero
	m_J.resize(numDoFs, numDoFs);
	m_J.set(0.0);

//	assemble matrix (depending on u, i.e. J(u))
	if(m_pAss->assemble_jacobian(m_J, u, *m_pDoFDistribution) != IAssemble_OK)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::init:"
				" Cannot assemble Jacobi matrix.\n");
		return false;
	}

//	Remember parallel storage type
	#ifdef UG_PARALLEL
		m_J.set_storage_type(PST_ADDITIVE);
		IDoFDistribution<TDoFDistribution>* dist =
				const_cast<IDoFDistribution<TDoFDistribution>*>(m_pDoFDistribution);
		CopyLayoutsAndCommunicatorIntoMatrix(m_J, *dist);
	#endif

//	remember that operator is initialized
	m_bInit = true;

//  we're done
	return true;
}

//	Initialize the operator
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
init()
{
//	todo: check that assembling is linear

//	check if DoF Distribution is set
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::prepare:"
				" DoF Distribution not set.\n");
		return false;
	}

//	get number of dofs
	const size_t numDoFs = m_pDoFDistribution->num_dofs();

//	Resize Matrix and set to zero
	m_J.resize(numDoFs, numDoFs);
	m_J.set(0.0);

//	Resize rhs
	m_rhs.resize(numDoFs);

//	Compute matrix (and rhs if needed)
	if(m_bAssembleRhs)
	{
	//	reset rhs to zero
		m_rhs.set(0.0);

	//	assemble matrix and rhs in one loop
		if(m_pAss->assemble_linear(m_J, m_rhs, m_rhs, *m_pDoFDistribution) != IAssemble_OK)
		{
			UG_LOG("ERROR in AssembledLinearOperator::init:"
					" Cannot assemble Matrix and Rhs.\n");
			return false;
		}

	//	remember parallel storage type of rhs
		#ifdef UG_PARALLEL
		m_rhs.set_storage_type(PST_ADDITIVE);
		#endif
	}
	else
	{
	//	assemble only matrix
		if(m_pAss->assemble_jacobian(m_J, m_rhs, *m_pDoFDistribution) != IAssemble_OK)
		{
			UG_LOG("ERROR in AssembledLinearOperator::init:"
					" Cannot assemble Matrix.\n");
			return false;
		}
	}

//	remember storage type for matrix
	#ifdef UG_PARALLEL
	m_J.set_storage_type(PST_ADDITIVE);
	IDoFDistribution<TDoFDistribution>* dist =
			const_cast<IDoFDistribution<TDoFDistribution>*>(m_pDoFDistribution);
	CopyLayoutsAndCommunicatorIntoMatrix(m_J, *dist);
	#endif

//	Remember that operator has been initialized
	m_bInit = true;

//	were done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
apply(vector_type& d, const vector_type& c)
{
	if(!m_bInit)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::apply: "
				"Operator not initialized.\n");
		return false;
	}

	#ifdef UG_PARALLEL
	if(!c.has_storage_type(PST_CONSISTENT))
		{
			UG_LOG("WARNING: In 'AssembledLinearizedOperator::apply':"
					"Inadequate storage format of Vector c.");
			return false;
		}
	#endif

	UG_ASSERT(c.size() == m_J.num_rows(),
	          "Row size '" << m_J.num_rows() << "' of Matrix J and size '"
	          << c.size() << "' of Vector x do not match. Cannot calculate L*x.");
	UG_ASSERT(d.size() == m_J.num_cols(),
	          "Column size '" << m_J.num_rows() << "' of Matrix J and size  '"
	          << d.size() << "' of Vector b do not match. Cannot calculate b := L*x.");

//	Apply Matrix
	return m_J.apply(d, c);
}

//	Compute d := d - J(u)*c
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
apply_sub(vector_type& d, const vector_type& c)
{
	if(!m_bInit)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::apply_sub: "
				"Operator not initialized.\n");
		return false;
	}

#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'AssembledLinearizedOperator::apply_sub':"
				"Inadequate storage format of Vector d.\n");
		return false;
	}
	if(!c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'AssembledLinearizedOperator::apply_sub':"
				"Inadequate storage format of Vector c.\n");
		return false;
	}
#endif

	UG_ASSERT(c.size() == m_J.num_rows(),
	          "Row size '" << m_J.num_rows() << "' of Matrix J and size '"
	          << c.size() << "' of Vector x do not match. Cannot calculate L*x.");
	UG_ASSERT(d.size() == m_J.num_cols(),
	          "Column size '" << m_J.num_rows() << "' of Matrix J and size  '"
	          << d.size() << "' of Vector b do not match. Cannot calculate b := b - L*x.");

//	Apply Matrix
	return m_J.matmul_minus(d,c);
}


template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
set_dirichlet_values(vector_type& u)
{
	if(!m_bInit)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::set_dirichlet_values:"
				" Operator not initialized.\n");
		return false;
	}

	if(m_pAss->assemble_solution(u, *m_pDoFDistribution) != IAssemble_OK)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" Cannot assemble solution.\n");
		return false;
	}
	return true;
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__ */

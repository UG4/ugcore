/*
 * assembled_linear_operator_impl.h
 *
 *  Created on: 24.02.2011
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__
#define __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__

#include "assembled_linear_operator.h"

//#define PROFILE_LIN_ASS
#ifdef PROFILE_LIN_ASS
	#define LINASS_PROFILE_FUNC()		PROFILE_FUNC()
	#define LINASS_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define LINASS_PROFILE_END()		PROFILE_END()
#else
	#define LINASS_PROFILE_FUNC()
	#define LINASS_PROFILE_BEGIN(name)
	#define LINASS_PROFILE_END()
#endif


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

//	assemble matrix (depending on u, i.e. J(u))
	if(m_pAss->assemble_jacobian(m_J, u, *m_pDoFDistribution) != true)
	{
		UG_LOG("ERROR in AssembledLinearizedOperator::init:"
				" Cannot assemble Jacobi matrix.\n");
		return false;
	}

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

//	resize rhs, since used as u dummy
	m_rhs.resize(m_pDoFDistribution->num_dofs());

//	Compute matrix (and rhs if needed)
	if(m_bAssembleRhs)
	{
	//	assemble matrix and rhs in one loop
		if(m_pAss->assemble_linear(m_J, m_rhs, m_rhs, *m_pDoFDistribution) != true)
		{
			UG_LOG("ERROR in AssembledLinearOperator::init:"
					" Cannot assemble Matrix and Rhs.\n");
			return false;
		}

	}
	else
	{
	//	assemble only matrix
		if(m_pAss->assemble_jacobian(m_J, m_rhs, *m_pDoFDistribution) != true)
		{
			UG_LOG("ERROR in AssembledLinearOperator::init:"
					" Cannot assemble Matrix.\n");
			return false;
		}
	}

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

	if(m_pAss->assemble_solution(u, *m_pDoFDistribution) != true)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" Cannot assemble solution.\n");
		return false;
	}
	return true;
}


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__OPERATOR__LINEAR_OPERATOR__ASSEMBLED_LINEAR_OPERATOR_IMPL__ */

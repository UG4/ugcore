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
		UG_LOG("ERROR in AssembledLinearOperator::init: "
				"DoF Distribution not set.\n");
		return false;
	}

	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Assembling routine not set.\n");
		return false;
	}

//	assemble matrix (depending on u, i.e. J(u))
	if(!m_pAss->assemble_jacobian(m_J, u, *m_pDoFDistribution))
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Cannot assemble Jacobi matrix.\n");
		return false;
	}

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
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" DoF Distribution not set.\n");
		return false;
	}
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::init:"
				" Assembling routine not set.\n");
		return false;
	}

//	resize rhs, since used as u dummy
	m_rhs.resize(m_pDoFDistribution->num_dofs());

//	Compute matrix (and rhs if needed)
	if(m_bAssembleRhs)
	{
	//	assemble matrix and rhs in one loop
		if(!m_pAss->assemble_linear(m_J, m_rhs, m_rhs, *m_pDoFDistribution))
		{
			UG_LOG("ERROR in AssembledLinearOperator::init:"
					" Cannot assemble Matrix and Rhs.\n");
			return false;
		}

	}
	else
	{
	//	assemble only matrix
		if(!m_pAss->assemble_jacobian(m_J, m_rhs, *m_pDoFDistribution))
		{
			UG_LOG("ERROR in AssembledLinearOperator::init:"
					" Cannot assemble Matrix.\n");
			return false;
		}
	}

//	were done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
apply(vector_type& d, const vector_type& c)
{
#ifdef UG_PARALLEL
	if(!c.has_storage_type(PST_CONSISTENT))
		{
			UG_LOG("WARNING: In 'AssembledLinearOperator::apply':"
					"Inadequate storage format of Vector c.");
			return false;
		}
#endif

	if(c.size() != m_J.num_cols() || d.size() == m_J.num_rows())
	{
		UG_LOG("ERROR in 'AssembledLinearOperator::apply': Size of matrix A ["<<
		        m_J.num_rows() << " x " << m_J.num_cols() << "] must match the "
		        "sizes of vectors x ["<<c.size()<<"], b ["<<d.size()<<"] for the "
		        " operation b = A*x. Maybe the operator is not initialized ?\n";)
		return false;
	}

//	Apply Matrix
	return m_J.apply(d, c);
}

//	Compute d := d - J(u)*c
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
apply_sub(vector_type& d, const vector_type& c)
{
#ifdef UG_PARALLEL
	if(!d.has_storage_type(PST_ADDITIVE))
	{
		UG_LOG("ERROR: In 'AssembledLinearOperator::apply_sub':"
				"Inadequate storage format of Vector d.\n");
		return false;
	}
	if(!c.has_storage_type(PST_CONSISTENT))
	{
		UG_LOG("ERROR: In 'AssembledLinearOperator::apply_sub':"
				"Inadequate storage format of Vector c.\n");
		return false;
	}
#endif

//	check sizes
	if(c.size() != m_J.num_cols() || d.size() == m_J.num_rows())
	{
		UG_LOG("ERROR in AssembledLinearOperator::apply_sub: Size of matrix A ["<<
		        m_J.num_rows() << " x " << m_J.num_cols() << "] must match the "
		        "sizes of vectors x ["<<c.size()<<"], b ["<<d.size()<<"] for the "
		        " operation b -= A*x. Maybe the operator is not initialized ?\n";)
		return false;
	}

//	Apply Matrix
	return m_J.matmul_minus(d,c);
}


template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledLinearOperator<TDoFDistribution, TAlgebra>::
set_dirichlet_values(vector_type& u)
{
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" Assembling routine not set.\n");
		return false;
	}
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in AssembledLinearOperator::set_dirichlet_values:"
				" DoF Distribution not set.\n");
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

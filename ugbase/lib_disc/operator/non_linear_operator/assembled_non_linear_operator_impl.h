/*
 * assembled_non_linear_operator_impl.h
 *
 *  Created on: ..
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__

#include "assembled_non_linear_operator.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledOperator<TDoFDistribution, TAlgebra>::
init()
{
	if(m_pDoFDistribution == NULL)
	{
		UG_LOG("ERROR in 'AssembledOperator::init':"
				" DoF Distribution not set.\n");
		return false;
	}
	if(m_pAss == NULL)
	{
		UG_LOG("ERROR in 'AssembledOperator::init': "
				"Discretization not set.\n");
		return false;
	}

//	remember that operator has been init
	m_bInit = true;

	return true;
}

//	Prepare functions
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledOperator<TDoFDistribution, TAlgebra>::
prepare(vector_type& dOut, vector_type& uIn)
{
	if(!m_bInit)
	{
		UG_LOG("ERROR in 'AssembledOperator::prepare': "
				"Operator not initialized.\n");
		return false;
	}

// 	Set Dirichlet - Nodes to exact values
	if(m_pAss->assemble_solution(uIn, *m_pDoFDistribution) != true)
	{
		UG_LOG("ERROR in 'AssembledOperator::prepare': "
				"Cannot set dirichlet values in solution.\n");
		return false;
	}

//	we're done
	return true;
}

// 	Compute d = L(u)
template <typename TDoFDistribution, typename TAlgebra>
bool
AssembledOperator<TDoFDistribution, TAlgebra>::
apply(vector_type& dOut, const vector_type& uIn)
{
	if(!m_bInit)
	{
		UG_LOG("ERROR in 'AssembledOperator::apply':"
				" Operator not initialized.\n");
		return false;
	}

//  assemble defect
	if(m_pAss->assemble_defect(dOut, uIn, *m_pDoFDistribution) != true)
	{
		UG_LOG("ERROR in 'AssembledOperator::apply': Could not "
				"assemble defect. Aborting.\n");
		return false;
	}

//	we're done
	return true;
}

} // end namepace ug

#endif /*__H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__*/

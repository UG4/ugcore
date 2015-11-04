
#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__

#include "assembled_non_linear_operator.h"

namespace ug{

//	Prepare functions
template <typename TAlgebra>
void
AssembledOperator<TAlgebra>::prepare(vector_type& uIn)
{
	PROFILE_BEGIN_GROUP(AssembledOperator_prepare, "discretization");
	if(m_spAss.invalid())
		UG_THROW("Discretization not set.");

// 	Set Dirichlet - Nodes to exact values (any constraints in general)
	try{
		m_spAss->adjust_solution(uIn, m_gridLevel);
	}
	UG_CATCH_THROW("Cannot set constraints in solution.");
}

// 	Compute d = L(u)
template <typename TAlgebra>
void
AssembledOperator<TAlgebra>::apply(vector_type& dOut, const vector_type& uIn)
{
	PROFILE_BEGIN_GROUP(AssembledOperator_apply, "discretization");
	if(m_spAss.invalid())
		UG_THROW("Discretization not set.");

//  assemble defect
	try{
		m_spAss->assemble_defect(dOut, uIn, m_gridLevel);
	}
	UG_CATCH_THROW("Could not assemble defect. Aborting.");
}

} // end namepace ug

#endif /*__H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__ASSEMBLED_NON_LINEAR_OPERATOR_IMPL__*/

#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__

#include "operator_interface.h"
#include "operator_inverse_interface.h"

namespace ug{

template <typename vector_type>
bool ApplyLinearSolver(	ILinearOperator<vector_type, vector_type>& A,
						vector_type& u, vector_type& b,
						ILinearOperatorInverse<vector_type, vector_type>& solver)
{
// step 1: Prepare: Assemble matrix
	PROFILE_BEGIN(assembleLinearMatrix);
	if(!A.init())
		{UG_LOG("ApplyLinearSolver: Cannot init Operator.\n"); return false;}
	PROFILE_END();

// step 2: Init Linear Inverse Operator
	PROFILE_BEGIN(initLinearSolver);
	if(!solver.init(A))
		{UG_LOG("ApplyLinearSolver: Cannot init Inverse operator.\n"); return false;}
	PROFILE_END();

// step 4: Apply Operator
	PROFILE_BEGIN(applyLinearSolver);
	if(!solver.apply_return_defect(u,b))
		{UG_LOG("ApplyLinearSolver: Cannot apply Inverse operator.\n"); return false;}
	PROFILE_END();

	return true;
}


}
#endif // __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__

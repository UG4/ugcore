#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__

#include "operator_interface.h"
#include "operator_inverse_interface.h"
#include "common/profiler/profiler.h"

namespace ug{

template <typename vector_type>
bool ApplyLinearSolver(	ILinearOperator<vector_type, vector_type>& A,
						vector_type& u, vector_type& b,
						ILinearOperatorInverse<vector_type, vector_type>& solver)
{
// step 1: Init Linear Inverse Operator
	PROFILE_BEGIN(ALS_InitLinearSolver);
	if(!solver.init(A))
	{
		UG_LOG("ApplyLinearSolver: Cannot init Inverse operator.\n");
		return false;
	}
	PROFILE_END();

// step 2: Apply Operator
	PROFILE_BEGIN(ALS_ApplyLinearSolver);
	if(!solver.apply_return_defect(u,b))
	{
		UG_LOG("ApplyLinearSolver: Cannot apply Inverse operator.\n");
		return false;
	}
	PROFILE_END();

//	done
	return true;
}


}
#endif // __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__

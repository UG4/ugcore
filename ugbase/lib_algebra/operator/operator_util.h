#ifndef __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__
#define __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__


#include "common/profiler/profiler.h"
#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/interface/linear_operator.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"

namespace ug{

template <typename vector_type>
bool ApplyLinearSolver(	SmartPtr<ILinearOperator<vector_type> > A,
						vector_type& u, vector_type& b,
						SmartPtr<ILinearOperatorInverse<vector_type> > solver)
{
// step 1: Init Linear Inverse Operator
	PROFILE_BEGIN_GROUP(ALS_InitLinearSolver, "algebra");
	if(!solver->init(A))
	{
		UG_LOG("ApplyLinearSolver: Cannot init Inverse operator.\n");
		return false;
	}
	PROFILE_END_(ALS_InitLinearSolver);

// step 2: Apply Operator
	PROFILE_BEGIN(ALS_ApplyLinearSolver);
	if(!solver->apply_return_defect(u,b))
	{
		UG_LOG("ApplyLinearSolver: Cannot apply Inverse operator.\n");
		return false;
	}
	PROFILE_END_(ALS_ApplyLinearSolver);

//	done
	return true;
}


}
#endif // __H__LIB_ALGEBRA__OPERATOR__OPERATOR_UTIL__

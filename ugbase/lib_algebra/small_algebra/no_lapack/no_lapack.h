#include "lu_decomp.h"
#include "nolapack_invert.h"

#ifndef __H__UG__CPU_ALGEBRA__NO_LAPACK_H__
#define __H__UG__CPU_ALGEBRA__NO_LAPACK_H__


template<typename A_type, typename B_type>
int GeneralizedEigenvalueProblem(A_type &A, A_type &X,
		B_type &lambda, A_type &B, bool bSortEigenvalues=false)
{
	UG_ASSERT(0, "GeneralizedEigenvalueProblem is only implemented for LAPACK at the moment");
	return 0;
}


#endif

#include "lu_decomp.h"
#include "../small_algebra.h"

#ifndef __H__UG__CPU_ALGEBRA__NO_LAPACK_H__
#define __H__UG__CPU_ALGEBRA__NO_LAPACK_H__


template<typename A_type, typename B_type>
int GeneralizedEigenvalueProblem(DenseMatrix<A_type> &A, DenseMatrix<A_type> &X,
		DenseVector<B_type> &lambda, DenseMatrix<A_type> &B, bool bSortEigenvalues=false)
{
	UG_ASSERT(0, "GeneralizedEigenvalueProblem is only implemented for LAPACK at the moment");
}


#endif

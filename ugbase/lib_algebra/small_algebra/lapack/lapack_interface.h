/**
 * \file lapack_interface.h
 *
 * \author Martin Rupp
 *
 * \date 10.08.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#ifndef __H__UG__MARTIN_ALGEBRA__LAPACK_INTERFACE_H__
#define __H__UG__MARTIN_ALGEBRA__LAPACK_INTERFACE_H__

#include "lapack.h"
#include "../small_algebra.h"

#include "lapack_invert.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

template<typename A_type, typename B_type>
int GeneralizedEigenvalueProblem(DenseMatrix<A_type> &A, DenseMatrix<A_type> &X,
		DenseVector<B_type> &lambda, DenseMatrix<A_type> &B, bool bSortEigenvalues=false)
{
	size_t N = A.num_rows();
	UG_ASSERT(N == A.num_cols() && N == B.num_rows() && N == B.num_cols(), "");

	B_type alphar; alphar.resize(N);
	B_type alphai; alphai.resize(N);
	B_type beta; beta.resize(N);


	int worksize = -1;
	double dWorksize = 0;

	int info;

	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, NULL, 0, &dWorksize, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], NULL, N, &X(0,0), N, &dWorksize, worksize);
	UG_ASSERT(info == 0, "gegv: failed to detect worksize");

	worksize = dWorksize;
	double *dwork = new double[worksize];
	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, NULL, 0, dwork, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], NULL, N, &X(0,0), N, dwork, worksize);
	UG_ASSERT(info == 0, "gegv: failed calculate");

	delete[] dwork;

	for(int i=0; i<N; i++)
	{
		lambda[i] = alphar[i]/beta[i];
		UG_ASSERT(dabs(alphai[i]) < 1e-10, "complex values");
	}

	if(bSortEigenvalues)
	{

		// bubblesort
		for(int i=N-1; i>0; i--)
		{
			for(int j=0; j<i; j++)
				if(lambda[j] > lambda[j+1])
				{
					swap(lambda[j], lambda[j+1]);
					for(size_t r = 0; r<X.num_rows(); r++)
						swap(X(r, j), X(r, j+1));
				}
		}
	}

	return info;
}





}
#endif

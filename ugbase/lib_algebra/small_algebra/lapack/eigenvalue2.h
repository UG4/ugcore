/**
 * \file eigenvalue2.h
 *
 * \author Martin Rupp
 *
 * \date 27.03.2013
 *
 * Goethe-Center for Scientific Computing 2013.
 *
 *
 *
 */

#ifndef __H__UG__EIGENVALUE__
#define __H__UG__EIGENVALUE__

#ifdef LAPACK_AVAILABLE

#include "lapack.h"
#include "../small_algebra.h"

#include "lapack_invert.h"
#include "lapack_densematrix_inverse.h"
#include "lapack_interface.h"
#include <complex>
#include "../small_matrix/densematrix.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// GeneralizedEigenvalueProblem
//-------------------------------
/**
 *	solves the generalized eigenvalue problem
 *	\f$ A x = \lambda B x\f$ with the LAPACK function \sa gegv.
 *	Works for double, floats and complex values.
 * \param	A					DenseMatrix A
 * \param	X					DenseMatrix X of eigenvectors. Column X[., i] is eigenvector \f$ x_i: \; A x_i =\lambda_i b\f$.
 * \param	lambda				\f$\lambda\f$-DenseVector
 * \param	B					DenseMatrix B
 * \param	bSortEigenvalues	sorts the eigenvalues (experimental)
 *
 */
template<typename A_type, typename TLambdaVectorType>
int GeneralizedEigenvalueProblemComplex(DenseMatrix<A_type> &A, DenseMatrix<A_type> &X,
		TLambdaVectorType &lambda, DenseMatrix<A_type> &B, bool bSortEigenvalues=false)
{
	size_t N = A.num_rows();
	UG_ASSERT(N == A.num_cols() && N == B.num_rows() && N == B.num_cols(), "");

	std::vector<double> alphar; alphar.resize(N);
	std::vector<double> alphai; alphai.resize(N);
	std::vector<double> beta; beta.resize(N);


	int worksize = -1;
	double dWorksize = 0;

	int info;

	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, NULL, 0, &dWorksize, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], NULL, N, &X(0,0), N, &dWorksize, worksize);
	UG_ASSERT(info == 0, "gegv: failed to detect worksize");

	worksize = (int)dWorksize;
	double *dwork = new double[worksize];
	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, NULL, 0, dwork, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], NULL, N, &X(0,0), N, dwork, worksize);
	UG_ASSERT(info == 0, "gegv: failed calculate");

	delete[] dwork;

	for(size_t i=0; i<N; i++)
		lambda[i] = std::complex<double>(alphar[i]/beta[i], alphai[i]/beta[i]);

	if(bSortEigenvalues)
	{

		// bubblesort
		for(int i=N-1; i>0; i--)
		{
			for(int j=0; j<i; j++)
				if(std::abs(lambda[j]) > std::abs(lambda[j+1]))
				{
					std::swap(lambda[j], lambda[j+1]);
					for(size_t r = 0; r<X.num_rows(); r++)
						std::swap(X(r, j), X(r, j+1));
				}
		}
	}

	return info;
}


}
#else // LAPACK_AVAILABLE
namespace ug{
template<typename A_type, typename TLambdaVectorType>
int GeneralizedEigenvalueProblemComplex(DenseMatrix<A_type> &A, DenseMatrix<A_type> &X,
		TLambdaVectorType &lambda, DenseMatrix<A_type> &B, bool bSortEigenvalues=false)
{
	UG_ASSERT(0, "GeneralizedEigenvalueProblemComplex is only implemented for LAPACK at the moment");
	return 0;
}
}
#endif // LAPACK_AVAILABLE
#endif // __H__UG__EIGENVALUE__

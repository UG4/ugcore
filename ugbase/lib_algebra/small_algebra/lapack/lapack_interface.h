/*
 * Copyright (c) 2010-2014:  G-CSC, Goethe University Frankfurt
 * Author: Martin Rupp
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

/**
 * \file lapack_interface.h
 *
 * \author Martin Rupp
 *
 * \date 10.08.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 *
 *
 * This file provides sophisticated functions for using DenseMatrices with LAPACK.
 *
 */

#ifndef __H__UG__MARTIN_ALGEBRA__LAPACK_INTERFACE_H__
#define __H__UG__MARTIN_ALGEBRA__LAPACK_INTERFACE_H__

#ifdef LAPACK_AVAILABLE

#include "lapack.h"
#include "../small_algebra.h"

#include "lapack_invert.h"
#include "lapack_densematrix_inverse.h"

namespace ug{

extern "C" // solve generalized eigenvalue system
	void dggev_(char *cCalcLeftEV, char *cCalcRightEV, lapack_int *n,
			lapack_double *pColMajorMatrixA, lapack_int *lda, lapack_double *pColMajorMatrixB,
			lapack_int *ldb, lapack_double *alphar,	lapack_double *alphai, lapack_double *beta,
			lapack_double *vl, lapack_int *ldvl, lapack_double *vr, lapack_int *ldvr, lapack_double *pWork,
				lapack_int *worksize, lapack_int *info);

inline lapack_int gegv(bool calcLeftEV, bool calcRightEV, lapack_int n,
		lapack_double *pColMajorMatrixA, lapack_int leading_dim_a, lapack_double *pColMajorMatrixB, lapack_int leading_dim_b,
		lapack_double *alphar, lapack_double *alphai, lapack_double *beta, lapack_double *vl, lapack_int leading_dim_vl,
		lapack_double *vr, lapack_int leading_dim_vr, lapack_double *work, lapack_int worksize)
{
	char jobvl = calcLeftEV ? 'V' : 'N';
	char jobvr = calcRightEV ? 'V' : 'N';
	lapack_int info=0;
	dggev_(&jobvl, &jobvr, &n, pColMajorMatrixA, &leading_dim_a, pColMajorMatrixB, &leading_dim_b, alphar, alphai, beta, vl, &leading_dim_vl, vr, &leading_dim_vr, work, &worksize, &info);
	return info;
}

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
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, nullptr, 0, &dWorksize, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], nullptr, N, &X(0,0), N, &dWorksize, worksize);
	UG_COND_THROW(info != 0, "gegv: failed to detect worksize");

	worksize = (int)dWorksize;
	double *dwork = new double[worksize];
	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, nullptr, 0, dwork, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], nullptr, N, &X(0,0), N, dwork, worksize);
	UG_COND_THROW(info != 0, "gegv: failed calculate");

	delete[] dwork;

	for(size_t i=0; i<N; i++)
	{
		lambda[i] = alphar[i]/beta[i];
		UG_COND_THROW(dabs(alphai[i]) > 1e-11, "gegv: complex values");
	}

	if(bSortEigenvalues)
	{

		// bubblesort
		for(int i=N-1; i>0; i--)
		{
			for(int j=0; j<i; j++)
				if(lambda[j] > lambda[j+1])
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

#endif // LAPACK_AVAILABLE
#endif

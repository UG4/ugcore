/*
 * Copyright (c) 2013-2014:  G-CSC, Goethe University Frankfurt
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
#include "common/error.h"

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
	THROW_IF_NOT_EQUAL_4(N, A.num_cols(), B.num_rows(), B.num_cols());

	std::vector<double> alphar; alphar.resize(N);
	std::vector<double> alphai; alphai.resize(N);
	std::vector<double> beta; beta.resize(N);


	int worksize = -1;
	double dWorksize = 0;

	int info;

	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, nullptr, 0, &dWorksize, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], nullptr, N, &X(0,0), N, &dWorksize, worksize);
	UG_COND_THROW(info != 0, "gegv: failed to detect worksize");

	worksize = static_cast<int>(dWorksize);
	auto *dwork = new double[worksize];
	if(A_type::ordering == RowMajor)
		info = gegv(true, false, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], &X(0,0), N, nullptr, 0, dwork, worksize);
	if(A_type::ordering == ColMajor)
		info = gegv(false, true, N, &A(0,0), N, &B(0,0), N, &alphar[0], &alphai[0], &beta[0], nullptr, N, &X(0,0), N, dwork, worksize);
	UG_COND_THROW(info != 0, "gegv: failed calculate");

	delete[] dwork;

	for(size_t i=0; i<N; i++)
		lambda[i] = std::complex<double>(alphar[i]/beta[i], alphai[i]/beta[i]);

	if(bSortEigenvalues)
	{

		// bubblesort
		for(int i=N-1; i>0; --i)
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
#endif
#endif

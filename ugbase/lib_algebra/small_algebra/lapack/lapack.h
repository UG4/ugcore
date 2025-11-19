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
 * \file lapack.h
 *
 * \author Martin Rupp
 *
 * \date 03.08.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 *
 *
 * General Lapack bridge functionality.
 * - overloaded functions for all lapack functions like sgetrf/dgetrf/cgetrf -> getrf,
 * - documentation of the lapack functions \sa getrf.
 * - the functions dont use pointers for all values (for example, m and n).
 * - enum eTransposeMode for transpose information.
 * - if you are missing a lapack function, this is the place to put it in
 *
 */

#ifndef __H__UG__CPU_ALGEBRA__LAPACK_H__
#define __H__UG__CPU_ALGEBRA__LAPACK_H__

#include "common/common.h"
#include <complex>

/*
#ifdef LAPACK_AVAILABLE
#include <clapack.h>
#endif
*/
namespace ug {

/*#ifndef __CLPK_integer
using lapack_int = __CLPK_integer;
using lapack_float = __CLPK_real;
using lapack_double = __CLPK_doublereal;
using lapack_ftnlen = __CLPK_ftnlen;
#else */
using lapack_int = int;
using lapack_float = float;
using lapack_double = double;
using lapack_ftnlen = int;

//#endif

enum eTransposeMode
{
	ModeNoTrans = 0,
	ModeTranspose = 1,
	ModeConjTranspose = 2
};

extern "C"
{
	// factor system *GETRF
    void sgetrf_(lapack_int *m, lapack_int *n, lapack_float *pColMajorMatrix, lapack_int *lda, lapack_int *ipiv, lapack_int *info);
    void dgetrf_(lapack_int *m, lapack_int *n, lapack_double *pColMajorMatrix, lapack_int *lda, lapack_int *ipiv, lapack_int *info);
    void cgetrf_(lapack_int *m, lapack_int *n, std::complex<lapack_float> *pColMajorMatrix, lapack_int *lda, lapack_int *ipiv, lapack_int *info);
    void zgetrf_(lapack_int *m, lapack_int *n, std::complex<lapack_double> *pColMajorMatrix, lapack_int *lda, lapack_int *ipiv, lapack_int *info);

	// solve system *GETRS
    void sgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const lapack_float *pColMajorMatrix,
    		lapack_int *lda, const lapack_int *ipiv, lapack_float *b, lapack_int *ldb, lapack_int *info);
    void dgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const lapack_double *pColMajorMatrix,
    		lapack_int *lda, const lapack_int *ipiv, lapack_double *b, lapack_int *ldb, lapack_int *info);
    void cgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const std::complex<lapack_float> *pColMajorMatrix,
    		lapack_int *lda, const lapack_int *ipiv, std::complex<lapack_float> *b, lapack_int *ldb, lapack_int *info);
	void zgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const std::complex<lapack_double> *pColMajorMatrix,
			lapack_int *lda, const lapack_int *ipiv, std::complex<lapack_double> *b, lapack_int *ldb, lapack_int *info);
	
	// invert system
	void sgetri_(lapack_int *n, lapack_float *pColMajorMatrix, lapack_int *lda, const lapack_int *ipiv, 
				 lapack_float *pWork, lapack_int *worksize, lapack_int *info);
	void dgetri_(lapack_int *n, lapack_double *pColMajorMatrix, lapack_int *lda, const lapack_int *ipiv,
				 lapack_double *pWork, lapack_int *worksize, lapack_int *info);	
}


inline char TransposeModeToChar(eTransposeMode t, bool isComplex)
{
	switch(t)
	{
	case ModeNoTrans:
		return 'N';
	case ModeTranspose:
		return 'T';
	case ModeConjTranspose:
		return (isComplex ? 'C' : 'T');
	default:
		UG_THROW("wrong tranpose mode");
		return 'N';
	}
}


// factor system
//--------------------

/**
 *
 *  getrf computes an LU factorization of a general M-by-N matrix A
 *  using partial pivoting with row interchanges.
 *
 *  The factorization has the form
 *     A = P * L * U
 *  where P is a permutation matrix, L is lower triangular with unit
 *  diagonal elements (lower trapezoidal if m > n), and U is upper
 *  triangular (upper trapezoidal if m < n).
 *
 *
 *  \param rows The number of rows of the matrix A.  M >= 0.
 *  \param cols The number of columns of the matrix A.  N >= 0.
 *  \param pColMajorMatrix dimension (LDA,N)
 *          On entry, the M-by-N matrix to be factored.
 *          On exit, the factors L and U from the factorization
 *          A = P*L*U; the unit diagonal elements of L are not stored.
 *
 *  \param lda The leading dimension of the array A.  LDA >= max(1,M).
 *  \param pPivot array, dimension (min(M,N))
 *          The pivot indices; for 1 <= i <= min(M,N), row i of the
 *          matrix was interchanged with row IPIV(i).
 *  \return = 0:  successful exit
 *          < 0:  if = -i, the i-th argument had an illegal value
 *          > 0:  if = i, U(i,i) is exactly zero. The factorization
 *                has been completed, but the factor U is exactly
 *                singular, and division by zero will occur if it is used
 *                to solve a system of equations.*
 */
inline lapack_int getrf(lapack_int rows, lapack_int cols, lapack_float *pColMajorMatrix, lapack_int lda, lapack_int *pPivot)
{
    lapack_int info;
    sgetrf_(&rows, &cols, pColMajorMatrix, &lda, pPivot, &info);
    return info;
}

inline lapack_int getrf(lapack_int rows, lapack_int cols, lapack_double *pColMajorMatrix, lapack_int lda, lapack_int *pPivot)
{
    lapack_int info;
    dgetrf_(&rows, &cols, pColMajorMatrix, &lda, pPivot, &info);
    return info;
}

inline lapack_int getrf(lapack_int rows, lapack_int cols, std::complex<lapack_float> *pColMajorMatrix, lapack_int lda, lapack_int *pPivot)
{
    lapack_int info;
    cgetrf_(&rows, &cols, pColMajorMatrix, &lda, pPivot, &info);
    return info;
}


inline lapack_int getrf(lapack_int rows, lapack_int cols, std::complex<lapack_double> *pColMajorMatrix, lapack_int lda, lapack_int *pPivot)
{
    lapack_int info;
    zgetrf_(&rows, &cols, pColMajorMatrix, &lda, pPivot, &info);
    return info;
}


// solve system
//---------------

/*
 *  getrs solves a system of linear equations
 *     A * X = B  or  A' * X = B
 *  with a general N-by-N matrix A using the LU factorization computed
 *  by getrf.
 *
 *  \param	transposeMode
 *			Specifies the form of the system of equations:
 *          = NoTranspose :  A * X = B
 *          = Transpose:  A**T * X = B
 *          = ConjugateTranspose :  A**H * X = B
 *
 *  \param	n order of the matrix A.  N >= 0.
 *
 *  \param	nrOfRHS
 *			The number of right hand sides, i.e., the number of columns
 *          of the matrix B.  nrOfRHS >= 0.
 *
 *  \param	pColMajorMatrix the Matrix A in column major ordering.
 *			float/double/std::complex<float> or std::complex<double> array, dimension (lda,n)
 *          The factors L and U from the factorization A = P*L*U
 *          as computed by getrf.
 *
 *  \param	lda The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  \param  pPivot  The pivot indices from DGETRF; for 1<=i<=N, row i of the
 *          matrix was interchanged with row pPivot[i].
 *
 *  \param  pRHS array, dimension (ldb,nrOfRHS) (col major)
 *          On entry, the right hand side matrix B.
 *          On exit, the solution matrix X.
 *
 *  \param	ldb The leading dimension of the array B.  LDB >= max(1,N).
 *
 *  \return 0:  successful exit. < 0: -i, the i-th argument had an illegal value
 *
 */
inline lapack_int getrs(eTransposeMode transposeMode, lapack_int n, lapack_int nrOfRHS, const float *pColMajorMatrix, lapack_int lda,
		  const lapack_int *pPivot, lapack_float *pRHS, lapack_int ldb)
{
    lapack_int info;
    char _trans = TransposeModeToChar(transposeMode, false);
    sgetrs_(&_trans, &n, &nrOfRHS, pColMajorMatrix, &lda, pPivot, pRHS, &ldb, &info);
    return info;
}

inline lapack_int getrs(eTransposeMode transposeMode, lapack_int n, lapack_int nrOfRHS, const double *pColMajorMatrix, lapack_int lda,
		  const lapack_int *pPivot, lapack_double *pRHS, lapack_int ldb)
{
    lapack_int info;
    char _trans = TransposeModeToChar(transposeMode, false);
	dgetrs_(&_trans, &n, &nrOfRHS, pColMajorMatrix, &lda, pPivot, pRHS, &ldb, &info);
    return info;
}

/*lapack_int getrs(eTransposeMode transposeMode, lapack_int n, lapack_int nrOfRHS, const std::complex<lapack_float> *pColMajorMatrix, lapack_int lda,
		  const lapack_int *pPivot, const std::complex<lapack_float> *pRHS, lapack_int ldb)
{
    lapack_int info;
    char _trans = TransposeModeToChar(transposeMode, true);
    cgetrs_(&_trans, &n, &nrOfRHS, pColMajorMatrix, &lda, pPivot, pRHS, &ldb, &info);
    return info;
}

lapack_int getrs(eTransposeMode transposeMode, lapack_int n, lapack_int nrOfRHS, const std::complex<lapack_double> *pColMajorMatrix, lapack_int lda,
		  const lapack_int *pPivot, const std::complex<lapack_double> *pRHS, lapack_int ldb)
{
    lapack_int info;
    char _trans = TransposeModeToChar(transposeMode, true);
	zgetrs_(&_trans, &n, &nrOfRHS, pColMajorMatrix, &lda, pPivot, pRHS, &ldb, &info);
    return info;
}*/
	
/*
 *  getri computes the inverse of a matrix using the LU factorization
 *  computed by getrf.
 *  This method inverts U and then computes inv(A) by solving the system
 *  inv(A) *L = inv(U) for inv(A).
 *
 *  \param	n (in)
 *          The order of the matrix A.  N >= 0.
 *
 *  \param	pColMajorMatrix (in/out) double array, dimension (LDA,N)
 *          On entry, the factors L and U from the factorization
 *          A = P *L *U as computed by getrf.
 *          On exit, if return = 0, the inverse of the original matrix A.
 *
 *  \param	lda (in)
 *          The leading dimension of the array A.  LDA >= max(1,N).
 *
 *  \param	pPivot (in) integer array, dimension n
 *          The pivot indices from getrf; for 0<=i<N, row i of the
 *          matrix was interchanged with row pPivot(i).
 *
 *  \param	pWork (workspace/output) double array, dimension max(1, worksize)
 *          On exit, if return=0, then work[0] returns the optimal worksize.
 *
 *  \param	worksize (in)
 *          The dimension of the array pWork.  worksize >= max(1, n)
 *          For optimal performance worksize >= n * nb, where nb is
 *          the optimal blocksize returned by ilaenv.
 *
 *          If worksize = -1, then a workspace query is assumed; the routine
 *          only calculates the optimal size of the pWork array, returns
 *          this value as the first entry of the pWork array, and no error
 *          message related to worksize is issued by XERBLA.
 *
 *  \return = 0:  successful exit
 *          < 0:  if return = -i, the i-th argument had an illegal value
 *          > 0:  if return = i, U(i,i) is exactly zero; the matrix is
 *                singular and its inverse could not be computed.
 */
inline lapack_int getri(lapack_int n, lapack_float *pColMajorMatrix, lapack_int lda, const int *pPivot,
				 lapack_float *pWork, lapack_int worksize)
{
	int info;
	sgetri_(&n, pColMajorMatrix, &lda, pPivot, pWork, &worksize, &info);
	return info;
}
	
inline lapack_int getri(lapack_int n, lapack_double *pColMajorMatrix, lapack_int lda, const int *pPivot,
	  lapack_double *pWork, lapack_int worksize)
{
	int info;
	dgetri_(&n, pColMajorMatrix, &lda, pPivot, pWork, &worksize, &info);
	return info;
}

}


#include "lapack_interface.h"

#endif



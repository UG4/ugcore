
/*
 *  stub_lapack.c
 *
 *
 *  Created by Martin Rupp on 02.08.10.
 *  Copyright 2010 . All rights reserved.
 *
 */

#ifndef __H__UG__CPU_ALGEBRA__LAPACK_LU_H__
#define __H__UG__CPU_ALGEBRA__LAPACK_LU_H__

#ifndef LAPACK_AVAILABLE

#include "lapack.h"

// factor system *GETRF
void sgetrf_(lapack_int *m, lapack_int *n, lapack_float *a, lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}
void dgetrf_(lapack_int *m, lapack_int *n, lapack_double *a, lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}
void cgetrf_(lapack_int *m, lapack_int *n, complex<lapack_float> *a, lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}
void zgetrf_(lapack_int *m, lapack_int *n, complex<lapack_double> *a, lapack_int *lda, lapack_int *ipiv, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}

// solve system *GETRS
void sgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const lapack_float *a,
		lapack_int *lda, const lapack_int *ipiv, lapack_float *b, lapack_int *ldb, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}
void dgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const lapack_double *a,
  		lapack_int *lda, const lapack_int *ipiv, lapack_double *b, lapack_int *ldb, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}
void cgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const complex<lapack_float> *a,
  		lapack_int *lda, const lapack_int *ipiv, lapack_double *b, lapack_int *ldb, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}
void zgetrs_(char *trans, lapack_int *n, lapack_int *nrhs, const complex<lapack_double> *a,
		lapack_int *lda, const lapack_int *ipiv, lapack_double *b, lapack_int *ldb, lapack_int *info)
{
	UG_ASSERT(0, "lapack not available");
}


#endif // LAPACK_AVAILABLE

#endif // __H__UG__CPU_ALGEBRA__LAPACK_H__

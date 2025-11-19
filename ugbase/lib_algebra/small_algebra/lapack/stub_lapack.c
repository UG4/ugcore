/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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


#endif

#endif
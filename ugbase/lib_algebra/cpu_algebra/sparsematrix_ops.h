/*
 * Copyright (c) 2026:  Goethe University Frankfurt
 * Author: Martin Rupp, Arne Naegel 
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
#pragma once

#include "sparsematrix.h"
#include "vector.h"

#include "../algebra_common/connection.h"
#include "../algebra_common/matrixrow.h"
#include "../algebra_common/sparsematrix_util.h"
#include "../common/operations_mat/operations_mat.h" // for MatMult, MatMultAdd, MatMultTransposedAdd


namespace ug
{

/**
 * \brief Sets all diagonal entries to a and all off-diagonal entries to zero.
 *
 * @param A Sparse matrix to modify.
 * @param a Value assigned to diagonal entries.
 */
template<typename T>
void SparseMatrixSet(SparseMatrix<T> &A, double a)
{
	PROFILE_SPMATRIX(SparseMatrix_set);
	A.check_fragmentation();

    const size_t num_rows = A.num_rows();
    
    #pragma omp parallel for
	for(size_t row = 0; row < num_rows; ++row)
	{
		auto itEnd = A.end_row(row);
		for(auto it = A.begin_row(row); it != itEnd; ++it)
		{   
            it.value() = (it.index() == row) ? a : 0.0;
		}
	}
}

/**
 * \brief Scales all entries of a sparse matrix by d.
 *
 * @param A Sparse matrix to modify.
 * @param d Scaling factor.
 */
template<typename T>
void SparseMatrixScale(SparseMatrix<T> &A, double d)
{
    const size_t num_rows = A.num_rows();
    #pragma omp parallel for
	for(size_t row = 0; row < num_rows; ++row)
	{
		auto itEnd = A.end_row(row);
		for(auto it = A.begin_row(row); it != itEnd; ++it)
			it.value() *= d;
	}
}

/**
 * \brief Calculates dest = alpha1*v1 + beta1*A*w1.
 *
 * @param A Sparse matrix A.
 * @param dest Destination vector.
 * @param alpha1 Scaling factor for v1.
 * @param v1 Vector used to initialize dest.
 * @param beta1 Scaling factor for A*w1.
 * @param w1 Input vector for the matrix-vector multiplication.
 */
template<typename T, typename vector_t>
void SparseMatrixAxpy(const SparseMatrix<T> &A, vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1)
{
	PROFILE_SPMATRIX(SparseMatrix_axpy);
	A.check_fragmentation();

	const size_t num_rows = A.num_rows();
	if(alpha1 == 0.0)
	{
		#pragma omp parallel for
		for(size_t i = 0; i < num_rows; ++i)
		{
			auto itEnd = A.end_row(i);
			auto it = A.begin_row(i);
			if(it == itEnd)
			{
				dest[i] = 0.0;
				continue;
			}

			MatMult(dest[i], beta1, it.value(), w1[it.index()]);
			for(++it; it != itEnd; ++it)
				MatMultAdd(dest[i], 1.0, dest[i], beta1, it.value(), w1[it.index()]);
		}
	}
	else if(&dest == &v1)
	{
		#pragma omp parallel for
		for(size_t i = 0; i < num_rows; ++i)
		{
			if(alpha1 != 1.0)
				dest[i] *= alpha1;

			auto itEnd = A.end_row(i);
			for(auto it = A.begin_row(i); it != itEnd; ++it)
				MatMultAdd(dest[i], 1.0, dest[i], beta1, it.value(), w1[it.index()]);
		}
	}
	else
	{
		#pragma omp parallel for
		for(size_t i = 0; i < num_rows; ++i)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			auto itEnd = A.end_row(i);
			for(auto it = A.begin_row(i); it != itEnd; ++it)
				MatMultAdd(dest[i], 1.0, dest[i], beta1, it.value(), w1[it.index()]);
		}
	}
}

/**
 * \brief Calculates dest = alpha1*v1 + beta1*A^T*w1.
 *
 * @param A Sparse matrix A.
 * @param dest Destination vector.
 * @param alpha1 Scaling factor for v1.
 * @param v1 Vector used to initialize dest.
 * @param beta1 Scaling factor for A^T*w1.
 * @param w1 Input vector for the transposed matrix-vector multiplication.
 */
template<typename T, typename vector_t>
void SparseMatrixAxpyTransposed(const SparseMatrix<T> &A, vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1)
{
	PROFILE_SPMATRIX(SparseMatrix_axpy_transposed);
	A.check_fragmentation();

	if(&dest == &v1)
	{
		if(alpha1 == 0.0)
			dest.set(0.0);
		else if(alpha1 != 1.0)
			dest *= alpha1;
	}
	else if(alpha1 == 0.0)
		dest.set(0.0);
	else
		VecScaleAssign(dest, alpha1, v1);

	#pragma omp parallel for
	for(size_t i = 0; i < A.num_rows(); ++i)
	{
		auto itEnd = A.end_row(i);
		for(auto it = A.begin_row(i); it != itEnd; ++it)
		{
			if(it.value() != 0.0)
				MatMultTransposedAdd(dest[it.index()], 1.0, dest[it.index()], beta1, it.value(), w1[i]);
		}
	}
}



} // namespace ug


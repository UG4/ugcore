/*
 *  sparsematrix_impl.hpp
 *
 *  Created by Martin Rupp on 04.11.09.
 *  Copyright 2009 G-CSC, University of Frankfurt. All rights reserved.
 *
 */
#ifndef __H__UG__CRS_ALGEBRA__SPARSEMATRIX_IMPL__
#define __H__UG__CRS_ALGEBRA__SPARSEMATRIX_IMPL__

// creation etc

#include <fstream>
#include <cstring>

#include "lib_algebra/common/operations_vec.h"
#include "common/profiler/profiler.h"
#include <vector>
#include <algorithm>



namespace ug{

template<typename T>
CRSSparseMatrix<T>::CRSSparseMatrix()
{
	bNeedsValues = true;
}

template<typename T>
bool CRSSparseMatrix<T>::resize(size_t newRows, size_t newCols)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_resize, "algebra CRSSparseMatrix");
	rowStart.clear(); rowStart.resize(newRows+1, -1);
	rowMax.clear(); rowMax.resize(newRows);
	rowEnd.clear(); rowEnd.resize(newRows, -1);
	m_numCols = newCols;
	nnz = 0;

	values.clear();
	if(bNeedsValues) values.resize(newRows);
	maxValues = 0;
	return true;
}


template<typename T>
bool CRSSparseMatrix<T>::set_as_transpose_of(const SparseMatrix<value_type> &B, double scale)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_set_as_transpose_of, "algebra CRSSparseMatrix");
	resize(B.num_cols(), B.num_rows());
	/*rowStart.resize(B.num_cols(), 0);
	rowMax.resize(B.num_cols(), 0);
	rowEnd.resize(B.num_cols(), 0);
	nnz = B.nnz;
	bNeedsValues = B.bNeedsValues;
	if(bNeedsValues) values.resize(nnz);
	cols.resize(nnz);
	maxValues = B.maxValues;
	fragmented = 0;
	m_numCols = B.num_rows();
	size_t r, c;
	for(r=0; r<num_rows(); r++)
		for(const_row_iterator it = B.begin_row(r); it != B.end_row(r); ++it)
			rowMax[it.index()]++;
	rowEnd[0] = rowMax[0];
	rowEnd[0] = 0;
	for(c=1; c<B.num_cols(); c++)
	{
		rowStart[c] = rowMax[c-1];
		rowMax[c] = rowStart[c]+rowMax[c];
		rowEnd[c] = rowStart[c];
	}*/

	for(size_t r=0; r<B.num_rows(); r++)
		for(const_row_iterator it = B.begin_row(r); it != B.end_row(r); ++it)
			operator()(it.index(), r) = scale*it.value();
	// todo: sort rows
	return true;
}

template<typename T>
template<typename vector_t>
inline void CRSSparseMatrix<T>::mat_mult_add_row(size_t row, typename vector_t::value_type &dest, double alpha, const vector_t &v) const
{
	for(const_row_iterator conn = begin_row(row); conn != end_row(row); ++conn)
		MatMultAdd(dest, 1.0, dest, alpha, conn.value(), v[conn.index()]);
}


// calculate dest = alpha1*v1 + beta1*A*w1 (A = this matrix)
template<typename T>
template<typename vector_t>
bool CRSSparseMatrix<T>::axpy(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1) const
{
	PROFILE_BEGIN_GROUP(SparseMatrix_axpy, "algebra CRSSparseMatrix");
	if(alpha1 == 0.0)
	{
		for(size_t i=0; i < num_rows(); i++)
		{
			const_row_iterator conn = begin_row(i);
			if(conn  == end_row(i)) continue;
			MatMult(dest[i], beta1, conn.value(), w1[conn.index()]);
			for(++conn; conn != end_row(i); ++conn)
				// res[i] += conn.value() * x[conn.index()];
				MatMultAdd(dest[i], 1.0, dest[i], beta1, conn.value(), w1[conn.index()]);
		}
	}
	else if(&dest == &v1)
	{
		if(alpha1 != 1.0) {
			for(size_t i=0; i < num_rows(); i++)
			{
				dest[i] *= alpha1;
				mat_mult_add_row(i, dest[i], beta1, w1);
			}
		}
		else
			for(size_t i=0; i < num_rows(); i++)
				mat_mult_add_row(i, dest[i], beta1, w1);

	}
	else
	{
		for(size_t i=0; i < num_rows(); i++)
		{
			VecScaleAssign(dest[i], alpha1, v1[i]);
			mat_mult_add_row(i, dest[i], beta1, w1);
		}
	}
	return true;
}

// calculate dest = alpha1*v1 + beta1*A^T*w1 (A = this matrix)
template<typename T>
template<typename vector_t>
bool CRSSparseMatrix<T>::axpy_transposed(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const vector_t &w1) const
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_axpy_transposed, "algebra CRSSparseMatrix");
	
	if(&dest == &v1) {
		if(alpha1 == 0.0)
			dest.set(0.0);
		else if(alpha1 != 1.0)
			dest *= alpha1;
	}
	else if(alpha1 == 0.0)
		dest.set(0.0);
	else
		VecScaleAssign(dest, alpha1, v1);

	for(size_t i=0; i<num_rows(); i++)
	{
		for(const_row_iterator conn = begin_row(i); conn != end_row(i); ++conn)
		{
			if(conn.value() != 0.0)
				// dest[conn.index()] += beta1 * conn.value() * w1[i];
				MatMultTransposedAdd(dest[conn.index()], 1.0, dest[conn.index()], beta1, conn.value(), w1[i]);
		}
	}
	return true;
}

template<typename T>
bool CRSSparseMatrix<T>::set(double a)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_set, "algebra CRSSparseMatrix");
	if(a == 0.0)
	{
		for(size_t row=0; row<num_rows(); row++)
			for(row_iterator it = begin_row(row); it != end_row(row); ++it)
				it.value() = 0.0;
	}
	else
	{
		for(size_t row=0; row<num_rows(); row++)
			for(row_iterator it = begin_row(row); it != end_row(row); ++it)
			{
				if(it.index() == row)
					it.value() = a;
				else
					it.value() = 0.0;
			}
	}
	return true;
}


template<typename T>
inline bool CRSSparseMatrix<T>::is_isolated(size_t i) const
{
	UG_ASSERT(i < num_rows() && i >= 0, *this << ": " << i << " out of bounds.");

	for(const_row_iterator it = begin_row(i); it != end_row(i); ++it)
		if(it.index() != i && it.value() != 0.0)
			return false;
	return true;
}


template<typename T>
void CRSSparseMatrix<T>::set_matrix_row(size_t row, connection *c, size_t nr)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_set_matrix_row, "algebra CRSSparseMatrix");
	for(size_t i=0; i<nr; i++)
		operator()(row, c[i].iIndex) = c[i].dValue;
}

template<typename T>
void CRSSparseMatrix<T>::add_matrix_row(size_t row, connection *c, size_t nr)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_add_matrix_row, "algebra CRSSparseMatrix");
	for(size_t i=0; i<nr; i++)
		operator()(row, c[i].iIndex) += c[i].dValue;
}


template<typename T>
bool CRSSparseMatrix<T>::set_as_copy_of(const CRSSparseMatrix<T> &B, double scale)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_set_as_copy_of, "algebra CRSSparseMatrix");
	resize(B.num_rows(), B.num_cols());
	for(size_t i=0; i < B.num_rows(); i++)
		for(const_row_iterator it = B.begin_row(i); it != B.end_row(i); ++it)
			operator()(i, it.index()) = scale*it.value();
}



template<typename T>
bool CRSSparseMatrix<T>::scale(double d)
{
	PROFILE_BEGIN_GROUP(CRSSparseMatrix_scale, "algebra CRSSparseMatrix");
	for(size_t i=0; i < num_rows(); i++)
		for(row_iterator it = begin_row(i); it != end_row(i); ++it)
			it.value() *= d;
	return true;
}







} // namespace ug

#endif


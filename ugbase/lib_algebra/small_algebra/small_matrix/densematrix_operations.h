/**
 * \file densematrix_operations.h
 *
 * \author Martin Rupp
 *
 * \date 29.10.2010
 *
 * Goethe-Center for Scientific Computing 2010.
 */

#include "densematrix.h"
#include "densevector.h"

#ifndef __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__
#define __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__

template<typename T>
struct matrix_algebra_type_traits<DenseMatrix<T> >
{
	static const matrix_algebra_type type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

//! calculates dest = beta1 * A1;
template<typename vector_t, typename matrix_t>
inline void MatMult(DenseVector<vector_t> &dest,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		MatMult(dest[r], A1(r,0), v[r]);
		for(size_type c = 1; c < v.size(); ++c)
			MatMultAdd(dest[r], 1.0, dest[r], A1(r,c), w1[c]);
	}
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultAdd(DenseVector<vector_t> &dest,
		const number &alpha1, const DenseVector<vector_t> &v1,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		VecScale(dest[r], alpha1, v1);
		for(size_type c = 0; c < v.size(); ++c)
			MatMultAdd(dest[r], 1.0, dest[r], A1(r,c), w1[c]);
	}
}




#endif // __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__

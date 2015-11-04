
#ifndef __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__
#define __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__

#include "densematrix.h"
#include "densevector.h"

#include "../../common/operations.h"

namespace ug{

/// \addtogroup small_algebra
/// \{

template<typename T>
struct matrix_algebra_type_traits;

template<typename T>
struct matrix_algebra_type_traits<DenseMatrix<T> >
{
	static const int type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

//! calculates dest = beta1 * A1 * w1;
template<typename vector_t, typename matrix_t>
inline void MatMult(DenseVector<vector_t> &dest,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		MatMult(dest[r], beta1, A1(r,0), w1[0]);
		for(size_t c = 1; c < w1.size(); ++c)
			MatMultAdd(dest[r], 1.0, dest[r], beta1, A1(r,c), w1[c]);
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
		VecScaleAssign(dest[r], alpha1, v1[r]);
		for(size_t c = 0; c < w1.size(); ++c)
			MatMultAdd(dest[r], 1.0, dest[r], beta1, A1(r,c), w1[c]);
	}
}

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(DenseVector<vector_t> &dest,
		const number &alpha1, const DenseVector<vector_t> &v1,
		const number &beta1, const DenseMatrix<matrix_t> &A1, const DenseVector<vector_t> &w1)
{
	for(size_t r = 0; r < dest.size(); ++r)
	{
		VecScaleAssign(dest[r], alpha1, v1[r]);
		for(size_t c = 0; c < w1.size(); ++c)
			MatMultTransposedAdd(dest[r], 1.0, dest[r], beta1, A1(c,r), w1[c]);
	}
}

// end group small_algebra
/// \}

}

#endif // __H__UG__COMMON__DENSEMATRIX_OPERATIONS_H__

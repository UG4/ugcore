/**
 * \file sparsematrix.h
 *
 * \author Martin Rupp
 *
 * \date 29.10.2012
 *
 * Goethe-Center for Scientific Computing 2012
 */


//////////////////// [NOTE] ///////////////////////
// this is currently only a demonstration to use different algebras, it maps to cpu_algebra

#ifndef __H__UG__CRS_ALGEBRA__SPARSEMATRIX__
#define __H__UG__CRS_ALGEBRA__SPARSEMATRIX__

#include "math.h"
#include "common/common.h"
#include "../algebra_common/sparsematrix_util.h"
#include "../common/operations_mat/operations_mat.h"
#include "../cpu_algebra/sparsematrix.h"

namespace ug{

/// \addtogroup crs_algebra
///	@{


/** CRSSparseMatrix
 *  \brief sparse matrix for big, variable sparse matrices.
 *
 * \param T blocktype
 */
template<typename TValueType> class CRSSparseMatrix : public SparseMatrix<TValueType>
{
	/// constructor for empty SparseMatrix
	CRSSparseMatrix() {}
	/// destructor
	virtual ~CRSSparseMatrix () {}
};


template<typename T>
struct matrix_algebra_type_traits<CRSSparseMatrix<T> >
{
	enum{
		type=MATRIX_USE_ROW_FUNCTIONS
	};
};

//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline void MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const CRSSparseMatrix<matrix_t> &A1, const vector_t &w1)
{
	A1.axpy_transposed(dest, alpha1, v1, beta1, w1);
}

// end group crs_algebra
/// \}

} // namespace ug

#endif

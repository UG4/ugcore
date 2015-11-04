/*
 * blas_mat_interface.h
 *
 *  Created on: 22.03.2011
 *      Author: mrupp
 */

#ifndef BLAS_MAT_INTERFACE_H_
#define BLAS_MAT_INTERFACE_H_



//! calculates dest = beta1 * A1 * w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultScale(vector_t &dest,	const number &beta1, const matrix_t &A1, const vector_t &w1);

//! calculates dest = A1 * w1;
template<typename vector_t, typename matrix_t>
inline bool MatMult(vector_t &dest,	const matrix_t &A1, const vector_t &w1);


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	);

//! calculates dest += A1*w1
template<typename vector_t, typename matrix_t>
inline bool MatMultAppend(vector_t &dest,
		const matrix_t &A1, const vector_t &w1	);


//! calculates dest = beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultScaleTransposed(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1	);

//! calculates dest = A1^T *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultScaleTransposed(vector_t &dest,
		const matrix_t &A1, const vector_t &w1	);


//! calculates dest = alpha1*v1 + beta1 * A1^T *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultTransposedScaledAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultTransposedAdd(dest, alpha1, v1, beta1, A1, w1);
}

// row functions

//! calculates dest += beta1 * A1[row] *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAddRow(typename vector_t::value_type &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1);
}

#endif /* BLAS_MAT_INTERFACE_H_ */

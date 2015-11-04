
#ifndef __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__
#define __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__

#include "common/common.h"
#include "matrix_algebra_types.h"

// operations for matrices
//-----------------------------------------------------------------------------
// these functions execute matrix operations

#include "matrix_use_operators.h"

#include "matrix_use_global_functions.h"

#include "matrix_use_row_functions.h"

#include "matrix_use_member_functions.h"

namespace ug
{

// functions transforming MatMult into the right operation, depending on matrix_algebra_type_traits.
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//! calculates dest = beta1 * A1;
template<typename vector_t, typename matrix_t>
inline bool MatMult(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
		::MatMult(dest, beta1, A1, w1);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1);
}

//! calculates dest = alpha1*v1 + alpha2*v2 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &alpha2, const vector_t &v2,
		const number &beta1, const matrix_t &A1, const vector_t &w1)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, alpha2, v2, beta1, A1, w1);
}

//! calculates dest = beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, beta1, A1, w1, A2, w2);
}

//! calculates dest = alpha1*v1 + beta1 * A1 *w1 + beta2 * A2*w2;
template<typename vector_t, typename matrix_t>
inline bool MatMultAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1,
		const number &beta2, const matrix_t &A2, const vector_t &w2)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultAdd(dest, alpha1, v1, beta1, A1, w1, beta2, A2, w2);
}


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultTransposed(vector_t &dest,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultTransposed(dest, beta1, A1, w1);
}


//! calculates dest = alpha1*v1 + beta1 * A1 *w1;
template<typename vector_t, typename matrix_t>
inline bool MatMultTransposedAdd(vector_t &dest,
		const number &alpha1, const vector_t &v1,
		const number &beta1, const matrix_t &A1, const vector_t &w1	)
{
	return mat_operations_class<vector_t, matrix_t, matrix_algebra_type_traits<matrix_t>::type>
			::MatMultTransposedAdd(dest, alpha1, v1, beta1, A1, w1);
}
/*
//! calculates dest = beta1*A1*w1 + alpha1*v1, and norm = norm^2(dest)
template<typename vector_t, typename matrix_t>
inline void MatMultAddNorm(vector_t &dest,
		const number &beta1, const TE_MAT_RowApplicable<matrix_t> &A1, const vector_t &w1,
		const number &alpha1, const vector_t &v1,
		double &norm)
{
	if(OPERATIONS_VERBOSE) cout << "dest: " << dest << endl <<
			"\t= " << beta1 << " * " << A1.cast() << "\t * " << w1 << endl <<
			"\t+ " << alpha1 << " * " << v1 << endl;

	norm=0;
	for(size_t i=0; i<dest.size(); i++)
	{
		VecScaleAssign(dest[i], alpha1, v1[i]);
		A1.cast().mat_mult_add_row(i, dest[i], beta1, w1);
		VecNormSquaredAdd(dest[i], norm);
	}
}
*/



} // namespace ug
#endif /* __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__ */

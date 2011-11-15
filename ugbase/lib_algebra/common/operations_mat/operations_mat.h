/*
 * operations_mat.h
 *
 *  Created on: 29.09.2010
 *      Author: mrupp
 */

#ifndef __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__
#define __H__UG__CPU_ALGEBRA__OPERATIONS_MAT__

#include "common/common.h"

namespace ug
{

// operations for matrices
//-----------------------------------------------------------------------------
// these functions execute matrix operations


/*
suppose you have your algebra with your vector class yourvector and your matrix class yourmatrix.
when your matrix class is rowwise multiplicable, and you want to do so, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static const matrix_algebra_type type = MATRIX_USE_ROW_FUNCTIONS;
};
and add a function
inline void yourmatrix::mat_mult_add_row(size_t row, yourvector::value_type &dest, double beta1, const yourvector &v) const;
to your class.

If you cannot or dont want to use rowwise multiplication, set
template<>
class matrix_algebra_type_traits<yourclass>
{
	static const matrix_algebra_type type = MATRIX_USE_GLOBAL_FUNCTIONS;
};

and implement (at least)
inline void MatMultAddDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1);
inline void MatMultDirect(yourvector &dest, const number &beta1, const yourmatrix &A1, const yourvector &w1,
							number &alpha1, const yourvector &v1);
Other functions like MatMultAdd(dest, beta1, A1, w1, beta2, A2, w2, alpha1, v1) are then implemented through your functions.
(see class mat_operations_class<vector_t, matrix_t, MATRIX_USE_GLOBAL_FUNCTIONS> for all functions).

*/


enum matrix_algebra_type
{
	MATRIX_USE_ROW_FUNCTIONS,
	MATRIX_USE_GLOBAL_FUNCTIONS,
	MATRIX_USE_OPERATORS,
	MATRIX_USE_MEMBER_FUNCTIONS
};

template<typename vector_t, typename matrix_t, int type>
struct mat_operations_class;

template<typename T>
struct matrix_algebra_type_traits
{
	static const int type = MATRIX_USE_OPERATORS;
};

#include "matrix_use_operators.h"

#include "matrix_use_global_functions.h"

#include "matrix_use_row_functions.h"

#include "matrix_use_member_functions.h"

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

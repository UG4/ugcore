/*
 * lgmath_matrix_vector_functions.h
 *
 *  Created on: 07.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__UGMATH__MATRIX_VECTOR_FUNCTIONS__
#define __H__UGMATH__MATRIX_VECTOR_FUNCTIONS__

#include "math_matrix.h"
#include "math_vector.h"

namespace ug{

/// Matrix - Vector Muliplication
// vOut = m * v
template <int N, int M>
inline
void
MatVecMult(MathVector<N>& vOut, const MathMatrix<N,M>& m, const MathVector<M>& v);

/// Transposed Matrix - Vector Muliplication
// vOut = Transpose(m) * v
template <int N, int M>
inline
void
TransposedMatVecMult(MathVector<N>& vOut, const MathMatrix<M,N>& m, const MathVector<M>& v);


} //end of namespace: lgmath

#include "math_matrix_vector_functions_common_impl.hpp"

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS__ */

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
template <typename vector_t, typename matrix_t>
inline
void
MatVecMult(vector_t& vOut, const matrix_t& m, const vector_t& v);

/// Transposed Matrix - Vector Muliplication
// vOut = Transpose(m) * v
template <typename vector_t, typename matrix_t>
inline
void
TransposedMatVecMult(vector_t& vOut, const matrix_t& m, const vector_t& v);


} //end of namespace: lgmath

#include "math_matrix_vector_functions_common_impl.hpp"

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS__ */

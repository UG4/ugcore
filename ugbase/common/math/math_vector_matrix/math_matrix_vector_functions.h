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

/// \addtogroup ugbase_math
/// \{

/// Matrix - Vector Multiplication
// vOut = m * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
MatVecMult(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v);

/// Matrix - Vector Multiplication adding to a second vector
// vOut += m * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
MatVecMultAppend(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v);

/// Matrix - Vector Multiplication added scaled to a second vector
// vOut += s * m * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
MatVecScaleMultAppend(vector_t_out& vOut, typename vector_t_out::value_type s, const matrix_t& m, const vector_t_in& v);

/// Transposed Matrix - Vector Muliplication
// vOut = Transpose(m) * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
TransposedMatVecMult(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v);

/// Transposed Matrix - Vector Muliplication
// vOut += Transpose(m) * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
TransposedMatVecMultAdd(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v);

// end group ugbase_math
/// \}

} //end of namespace: ug

#include "math_matrix_vector_functions_common_impl.hpp"

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS__ */

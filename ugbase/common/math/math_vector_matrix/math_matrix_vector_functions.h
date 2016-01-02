/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

/// Multiplication by the Givens rotation of the QR-Decomposition
template <typename matrix_t, typename vector_t>
inline
void
GivensMatVecMult (matrix_t& A, vector_t& v);

/// Multiplication by the inverse using the Givens rotations
template <typename matrix_t, typename vector_t>
inline
void
InvMatVecMult_byGivens (matrix_t& A, vector_t& v);

/// Orthogonal projection
template <typename matrix_t, typename vector_t>
inline
void
OrthogProjectVec (vector_t& v, const matrix_t& A);

// end group ugbase_math
/// \}

} //end of namespace: ug

#include "math_matrix_vector_functions_common_impl.hpp"

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS__ */

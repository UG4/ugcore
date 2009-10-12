/*
 * lgmath_matrix_functions.h
 *
 *  Created on: 07.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__COMMON__MATHMATRIX_FUNCTIONS__
#define __H__COMMON__MATHMATRIX_FUNCTIONS__

#include "math_matrix.h"

namespace ug{
////////////////////////////////////////////////////////////////
// Addition of Matrices

///	adds two matrices and stores the result in a third one
// mOut = m1 + m2
template <int N, int M>
inline
void
MatAdd(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m1, const MathMatrix<N,M>& m2);

////////////////////////////////////////////////////////////////
// Subtraction of Matrices

///	subtracts m2 from m1 and stores the result in a mOut
// mOut = m1 - m2
template <int N, int M>
inline
void
MatSubtract(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m1, const MathMatrix<N,M>& m2);

////////////////////////////////////////////////////////////////
// Scaling of Matrices

///	scales a MathMatrix<N,M>
// mOut = s * m
template <int N, int M>
inline
void
MatScale(MathMatrix<N,M>& mOut, typename MathMatrix<N,M>::value_type s, const MathMatrix<N,M>& m);

////////////////////////////////////////////////////////////////
// Determinant of Matrix

/// Determinant of a MathMatrix<N,M>
inline
MathMatrix<2,2>::value_type
Determinant(const MathMatrix<2,2>& m);
inline
MathMatrix<3,3>::value_type
Determinant(const MathMatrix<3,3>& m);

////////////////////////////////////////////////////////////////
// Transposed of Matrix

/// transpose a MathMatrix<N,M>
template <int N, int M>
inline
void
Transpose(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m);

/// transpose a MathMatrix<N,M>, override original MathMatrix<N,M>
template <int N, int M>
inline
void
Transpose(MathMatrix<N,M>& m);

////////////////////////////////////////////////////////////////
// Inverse of Matrix

/// Inverse of a MathMatrix<N,M>
inline
void
Inverse(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m);
inline
void
Inverse(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m);

/// Inverse of a MathMatrix<N,M>
inline
void
Inverse(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m, MathMatrix<2,2>::value_type& det);
inline
void
Inverse(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m, MathMatrix<3,3>::value_type& det);

////////////////////////////////////////////////////////////////
// Inverse Transposed of Matrix

/// Transposed-Inverse of a MathMatrix<N,M> (= Inverse-Transposed of a MathMatrix<N,M>)
inline
void
InverseTransposed(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m);
inline
void
InverseTransposed(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m);

/// Transposed-Inverse of a MathMatrix<N,M> (= Inverse-Transposed of a MathMatrix<N,M>)
inline
void
InverseTransposed(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m, MathMatrix<2,2>::value_type& det);
inline
void
InverseTransposed(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m, MathMatrix<3,3>::value_type& det);

////////////////////////////////////////////////////////////////
// Scalar operations for Matrices

/// Set each matrix entry to a scalar (componentwise)
template <int N, int M>
inline
void
MatSet(MathMatrix<N,M>& mInOut, typename MathMatrix<N,M>::value_type s);

/// Set each diagonal of a matrix to a scalar (componentwise)
template <int N, int M>
inline
void
MatDiagSet(MathMatrix<N,M>& mInOut, typename MathMatrix<N,M>::value_type s);

/// Add a scalar to a vector (componentwise)
template <int N, int M>
inline
void
MatAdd(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s);

/// Subtract a scalar from a vector (componentwise)
template <int N, int M>
inline
void
MatSubtract(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s);

/// Devide a vector by a scalar (componentwise)
template <int N, int M>
inline
void
MatDevide(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s);

/// Multiply a vector by a scalar (componentwise)
template <int N, int M>
inline
void
MatMultiply(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s);

////////////////////////////////////////////////////////////////
// Norms for Matrices

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatFrobeniusNormSq(MathMatrix<N,M>& m);

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatFrobeniusNorm(MathMatrix<N,M>& m);

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatOneNorm(MathMatrix<N,M>& m);

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatInftyNorm(MathMatrix<N,M>& m);

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatMaxNorm(MathMatrix<N,M>& m);



} //end of namespace: lgmath

////////////////////////////////////////////////////////////////////////
//	include a general, but not very fast implementation of the declared methods above.
#include "math_matrix_functions_common_impl.hpp"

#endif /* __H__LGMATH__LGMATH_MathMatrix<N,M>_FUNCTIONS__ */

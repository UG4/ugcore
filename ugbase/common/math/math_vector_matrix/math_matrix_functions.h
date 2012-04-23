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
template <typename matrix_t>
inline
void
MatAdd(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2);

////////////////////////////////////////////////////////////////
// Subtraction of Matrices

///	subtracts m2 from m1 and stores the result in a mOut
// mOut = m1 - m2
template <typename matrix_t>
inline
void
MatSubtract(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2);

////////////////////////////////////////////////////////////////
// Multiplication of Matrices

///	multiply two matrices and stores the result in a third one
// mOut = m1 * m2
template <size_t N, size_t M, size_t L, typename T>
inline
void
MatMult(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1, const MathMatrix<L, M, T>& m2);

///	multiply three matrices and stores the result in a fourth one
// mOut = m1 * m2 * m3
template <size_t N, size_t M, size_t L, size_t P, typename T>
inline
void
MatMult(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1, const MathMatrix<L, P, T>& m2, const MathMatrix<P, M, T>& m3);

///	multiply two transposed matrices and stores the result in a third one
// mOut = m1^T * m2^T
template <size_t N, size_t M, size_t L, typename T>
inline
void
MatMultiplyTransposed(MathMatrix<N, M, T>& mOut, const MathMatrix<L, N, T>& m1, const MathMatrix<M, L, T>& m2);

///	multiply a transposed matrix with itself and stores the result in a second one
// mOut = m^T * m
template <size_t N, size_t M, typename T>
inline
void
MatMultiplyMTM(MathMatrix<N, N, T>& mOut, const MathMatrix<M, N, T>& m);

///	multiply a matrix with its transposed and stores the result in a second one
// mOut = m * m^T
template <size_t N, size_t M, typename T>
inline
void
MatMultiplyMMT(MathMatrix<M, M, T>& mOut, const MathMatrix<M, N, T>& m);

////////////////////////////////////////////////////////////////
// "Contraction" for Matrices (note: contraction is only known regarding tensors!)

template <typename matrix_t>
inline
typename matrix_t::value_type
MatContraction(const matrix_t& m1, const matrix_t& m2);

////////////////////////////////////////////////////////////////
// Scaling of Matrices

///	scales a matrix_t
// mOut = s * m
template <typename matrix_t>
inline
void
MatScale(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m);

///	scales a matrix_t and adds to result to a second matrix
// mOut += s * m
template <typename matrix_t>
inline
void
MatScaleAppend(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m);

////////////////////////////////////////////////////////////////
// Determinant of Matrix

/// Determinant of a matrix_t
template <typename T>
inline
typename MathMatrix<1,1,T>::value_type
Determinant(const MathMatrix<1,1,T>& m);

template <typename T>
inline
typename MathMatrix<2,2,T>::value_type
Determinant(const MathMatrix<2,2,T>& m);

template <typename T>
inline
typename MathMatrix<3,3,T>::value_type
Determinant(const MathMatrix<3,3,T>& m);

////////////////////////////////////////////////////////////////
// Transposed of Matrix

/// transpose a matrix_t
template <typename matrix_t>
inline
void
Transpose(matrix_t& mOut, const matrix_t& m);

/// transpose a matrix_t, override original matrix_t
template <typename matrix_t>
inline
void
Transpose(matrix_t& m);

////////////////////////////////////////////////////////////////
// Inverse of Matrix

/// Inverse of a matrix_t
template <typename T>
inline
void
Inverse(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m);
template <typename T>
inline
void
Inverse(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m);
template <typename T>
inline
void
Inverse(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m);

/// Inverse of a matrix_t
template <typename T>
inline
void
Inverse(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m, typename MathMatrix<1,1,T>::value_type& det);
template <typename T>
inline
void
Inverse(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m, typename MathMatrix<2,2,T>::value_type& det);
template <typename T>
inline
void
Inverse(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m, typename MathMatrix<3,3,T>::value_type& det);

////////////////////////////////////////////////////////////////
// Inverse Transposed of Matrix

/// Transposed-Inverse of a matrix_t (= Inverse-Transposed of a matrix_t)
template <typename T>
inline
void
InverseTransposed(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m);
template <typename T>
inline
void
InverseTransposed(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m);
template <typename T>
inline
void
InverseTransposed(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m);

/// Transposed-Inverse of a matrix_t (= Inverse-Transposed of a matrix_t)
template <typename T>
inline
void
InverseTransposed(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m, typename MathMatrix<1,1,T>::value_type& det);
template <typename T>
inline
void
InverseTransposed(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m, typename MathMatrix<2,2,T>::value_type& det);
template <typename T>
inline
void
InverseTransposed(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m, typename MathMatrix<3,3,T>::value_type& det);

////////////////////////////////////////////////////////////////
// Right-Inverse of Matrix
template <size_t N, size_t M, typename T>
inline
void
RightInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m);

////////////////////////////////////////////////////////////////
// Left-Inverse of Matrix
template <size_t N, size_t M, typename T>
inline
void
LeftInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m);

////////////////////////////////////////////////////////////////
// Trace of Matrix

/// Trace of a matrix_t
template <typename T>
inline
typename MathMatrix<1,1,T>::value_type
Trace(const MathMatrix<1,1,T>& m);

template <typename T>
inline
typename MathMatrix<2,2,T>::value_type
Trace(const MathMatrix<2,2,T>& m);

template <typename T>
inline
typename MathMatrix<3,3,T>::value_type
Trace(const MathMatrix<3,3,T>& m);

////////////////////////////////////////////////////////////////
// Scalar operations for Matrices

/// Set each matrix entry to a scalar (componentwise)
template <typename matrix_t>
inline
void
MatSet(matrix_t& mInOut, typename matrix_t::value_type s);

/// Set each diagonal of a matrix to a scalar (componentwise)
template <typename matrix_t>
inline
void
MatDiagSet(matrix_t& mInOut, typename matrix_t::value_type s);

/// Add a scalar to a matrix (componentwise)
template <typename matrix_t>
inline
void
MatAdd(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

/// Subtract a scalar from a matrix (componentwise)
template <typename matrix_t>
inline
void
MatSubtract(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

/// Devide a matrix by a scalar (componentwise)
template <typename matrix_t>
inline
void
MatDevide(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

/// Multiply a matrix by a scalar (componentwise)
template <typename matrix_t>
inline
void
MatMultiply(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

///	Fills the matrix with the identity matrix.
template <typename matrix_t>
inline
void
MatIdentity(matrix_t& mOut);

///	Fills the matrix with a matrix that rotates around the x-axis in 3 dimensions.
/**	matrix_t has to have at least 3 rows and 3 columns.
 *	Rotation is specified in radiants.*/
template <typename matrix_t>
inline
void
MatRotationX(matrix_t& mOut, typename matrix_t::value_type rads);

///	Fills the matrix with a matrix that rotates around the y-axis in 3 dimensions.
/**	matrix_t has to have at least 3 rows and 3 columns.
 *	Rotation is specified in radiants.*/
template <typename matrix_t>
inline
void
MatRotationY(matrix_t& mOut, typename matrix_t::value_type rads);

///	Fills the matrix with a matrix that rotates around the y-axis in 2 or 3 dimensions.
/**	matrix_t has to have at least 2 rows and 2 columns.
 *	Rotation is specified in radiants.*/
template <typename matrix_t>
inline
void
MatRotationZ(matrix_t& mOut, typename matrix_t::value_type rads);

////////////////////////////////////////////////////////////////
///	Creates a rotation matrix given yaw, pitch and roll in radiants.
template <typename matrix_t>
inline
void
MatRotationYawPitchRoll(matrix_t& mOut,
						typename matrix_t::value_type yaw,
						typename matrix_t::value_type pitch,
						typename matrix_t::value_type roll);
						
////////////////////////////////////////////////////////////////
// Norms for Matrices

template <typename matrix_t>
inline
typename matrix_t::value_type
MatFrobeniusNormSq(matrix_t& m);

template <typename matrix_t>
inline
typename matrix_t::value_type
MatFrobeniusNorm(matrix_t& m);

template <typename matrix_t>
inline
typename matrix_t::value_type
MatOneNorm(matrix_t& m);

template <typename matrix_t>
inline
typename matrix_t::value_type
MatInftyNorm(matrix_t& m);

template <typename matrix_t>
inline
typename matrix_t::value_type
MatMaxNorm(matrix_t& m);



} //end of namespace: lgmath

////////////////////////////////////////////////////////////////////////
//	include a general, but not very fast implementation of the declared methods above.
#include "math_matrix_functions_common_impl.hpp"

#endif /* __H__LGMATH__LGMATH_MATRIX_FUNCTIONS__ */

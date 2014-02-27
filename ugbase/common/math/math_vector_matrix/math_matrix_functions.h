/*
 * math_matrix_functions.h
 *
 *  Created on: 07.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__COMMON__MATH_MATRIX_FUNCTIONS__
#define __H__UG__COMMON__MATH_MATRIX_FUNCTIONS__

#include "math_matrix.h"

namespace ug{

/// \addtogroup math_matrix
/// \{

////////////////////////////////////////////////////////////////////////////////
// Addition of Matrices
////////////////////////////////////////////////////////////////////////////////

///	adds two matrices and stores the result in a third one
// mOut = m1 + m2
template <typename matrix_t>
inline void
MatAdd(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2);

////////////////////////////////////////////////////////////////////////////////
// Subtraction of Matrices
////////////////////////////////////////////////////////////////////////////////

///	subtracts m2 from m1 and stores the result in a mOut
// mOut = m1 - m2
template <typename matrix_t>
inline void
MatSubtract(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2);

////////////////////////////////////////////////////////////////////////////////
// Multiplication of Matrices
////////////////////////////////////////////////////////////////////////////////

///	multiply two matrices and stores the result in a third one
// mOut = m1 * m2
template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiply(MathMatrix<N, M, T>& mOut,
        const MathMatrix<N, L, T>& m1, const MathMatrix<L, M, T>& m2);

///	multiply three matrices and stores the result in a fourth one
// mOut = m1 * m2 * m3
template <size_t N, size_t M, size_t L, size_t P, typename T>
inline void
MatMultiply(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1,
        const MathMatrix<L, P, T>& m2, const MathMatrix<P, M, T>& m3);

///	multiply two transposed matrices and stores the result in a third one
// mOut = m1^T * m2^T
template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiplyTransposed(MathMatrix<N, M, T>& mOut,
                      const MathMatrix<L, N, T>& m1, const MathMatrix<M, L, T>& m2);

///	multiply a transposed matrix with itself and stores the result in a second one
// mOut = m^T * m
template <size_t N, size_t M, typename T>
inline void
MatMultiplyMTM(MathMatrix<N, N, T>& mOut, const MathMatrix<M, N, T>& m);

///	multiply a matrix with its transposed and stores the result in a second one
// mOut = m * m^T
template <size_t N, size_t M, typename T>
inline void
MatMultiplyMMT(MathMatrix<M, M, T>& mOut, const MathMatrix<M, N, T>& m);

///	multiply a matrix with the transposed of a second one and stores the result in mOut
// mOut = m1 * m2^T
template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiplyMBT(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1,
        const MathMatrix<M, L, T>& m2);

///	multiply the transposed of a matrix with a matrix and stores the result in mOut
// mOut = m1^T * m2
template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiplyMTB(MathMatrix<N, M, T>& mOut, const MathMatrix<L, N, T>& m1,
        const MathMatrix<L, M, T>& m2);

///	multiply a matrix m2 with a matrix m1 from left and its transpose m1^T from right
///	and stores the result in mOut
// mOut = m1 * m2 * m1^T
template <size_t N, size_t M, typename T>
inline void
MatMultiplyMBMT(MathMatrix<N, N, T>& mOut, const MathMatrix<N, M, T>& m1,
        const MathMatrix<M, M, T>& m2);

///	multiply a matrix m2 with a matrix m1 from right and its transpose m1^T from left
///	and stores the result in mOut
// mOut = m1^T * m2 * m1
template <size_t N, size_t M, typename T>
inline void
MatMultiplyMTBM(MathMatrix<N, N, T>& mOut, const MathMatrix<M, N, T>& m1,
        const MathMatrix<M, M, T>& m2);

////////////////////////////////////////////////////////////////////////////////
// "Contraction" for Matrices (note: contraction is only known regarding tensors!)
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline typename matrix_t::value_type
MatContraction(const matrix_t& m1, const matrix_t& m2);

////////////////////////////////////////////////////////////////////////////////
// Scaling of Matrices
////////////////////////////////////////////////////////////////////////////////

///	scales a matrix_t
// mOut = s * m
template <typename matrix_t>
inline void
MatScale(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m);

///	scales a matrix_t and adds to result to a second matrix
// mOut += s * m
template <typename matrix_t>
inline void
MatScaleAppend(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m);

////////////////////////////////////////////////////////////////////////////////
// Transposed of Matrix
////////////////////////////////////////////////////////////////////////////////

/// transpose a matrix
template <size_t N, size_t M, typename T>
inline void
Transpose(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m);

/// transpose a matrix_t, override original matrix_t
template <typename matrix_t>
inline void
Transpose(matrix_t& m);

////////////////////////////////////////////////////////////////////////////////
// Determinant of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Determinant of a matrix
/**
 * Returns the Determinate of a matrix.
 *
 * @param m 	Matrix
 * @return determinant of matrix
 */
/// \{
template <size_t N, typename T>
inline typename MathMatrix<N,N,T>::value_type
Determinant(const MathMatrix<N,N,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
Determinant(const MathMatrix<1,1,T>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
Determinant(const MathMatrix<2,2,T>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
Determinant(const MathMatrix<3,3,T>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Gram Determinant of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Gram Determinant of a matrix
/**
 * Returns Gram Determinant of a matrix, i.e. \f$ \det(M M^T) \f$.
 * @param m		Matrix
 * @return	gram determinant of matrix
 */
/// \{
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
GramDeterminant(const MathMatrix<N,M,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
GramDeterminant(const MathMatrix<1,1,T>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
GramDeterminant(const MathMatrix<2,2,T>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
GramDeterminant(const MathMatrix<3,3,T>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Square root of Gram Determinant of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Square root of Gram Determinant of a matrix
/**
 * Returns Square root of Gram Determinant of a matrix,
 * i.e. \f$ \sqrt{\det(M M^T)} \f$.
 * @param m		Matrix
 * @return	square root of gram determinant of matrix
 */
/// \{
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
SqrtGramDeterminant(const MathMatrix<N,M,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
SqrtGramDeterminant(const MathMatrix<1,1,T>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
SqrtGramDeterminant(const MathMatrix<2,2,T>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
SqrtGramDeterminant(const MathMatrix<3,3,T>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Inverse of a matrix
/**
 * Computes the inverse of a matrix and returns the determinante.
 *
 * @param mOut		Inverse of Matrix
 * @param m			Matrix
 * @return		Determinate of Matrix
 */
/// \{
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
Inverse(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
Inverse(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
Inverse(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
Inverse(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Inverse Transposed of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Transposed-Inverse of a Matrix (= Inverse-Transposed of a Matrix)
/**
 * Computes the Inverse-Transposed of a Matrix and returns the Determinant.
 *
 * @param mOut 		Inverse-Transposed of Matrix
 * @param m			Matrix
 * @return			Determinant of Matrix
 */
/// \{
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
InverseTransposed(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
InverseTransposed(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
InverseTransposed(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
InverseTransposed(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Right-Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Right-Inverse of a Matrix
/**
 * Computes the Right-Inverse of a Matrix and returns the square root of the
 * gram determinate. For any \f$ A \in \mathbb{R}^{M \times N} \f$
 *  (with \f$ M \leq N \f$) the Right-Inverse is defined as
 * \f[
 * 		A^{-1}_{right} := A^T (A A^T)^{-1},
 * \f]
 * such that \f$ A A^{-1}_{right} = 1 \f$. The gram determinate is defined as
 * \f$ \det{A A^T} \f$.
 *
 * Please note, that in case of a square matrix (i.e. \f$ N = M \f$), the Right-
 * Inverse is identical to the Inverse and the square root of the gram determinante
 * simplifies to the norm of the usual determinate.
 *
 * @param mOut 		Right-Inverse of Matrix
 * @param m			Matrix
 * @return			Square root of gram determinate of Matrix
 */
/// \{
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
RightInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
RightInverse(MathMatrix<1,1>& mOut, MathMatrix<1,1>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
RightInverse(MathMatrix<2,2>& mOut, MathMatrix<2,2>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
RightInverse(MathMatrix<3,3>& mOut, MathMatrix<3,3>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Left-Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Left-Inverse of a Matrix
/**
 * Computes the Left-Inverse of a Matrix and returns the square root of the
 * gram determinate. For any \f$ A \in \mathbb{R}^{M \times N} \f$
 *  (with \f$ M \geq N \f$) the Left-Inverse is defined as
 * \f[
 * 		A^{-1}_{left} := (A^T A)^{-1} A^T,
 * \f]
 * such that \f$ A^{-1}_{left} A = 1 \f$. The gram determinate is defined as
 * \f$ \det{A^T A} \f$.
 *
 * Please note, that in case of a square matrix (i.e. \f$ N = M \f$), the Left-
 * Inverse is identical to the Inverse and the square root of the gram determinante
 * simplifies to the norm of the usual determinate.
 *
 * @param mOut 		Left-Inverse of Matrix
 * @param m			Matrix
 * @return			Square root of gram determinate of Matrix
 */
/// \{
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
LeftInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m);

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
LeftInverse(MathMatrix<1,1>& mOut, MathMatrix<1,1>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
LeftInverse(MathMatrix<2,2>& mOut, MathMatrix<2,2>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
LeftInverse(MathMatrix<3,3>& mOut, MathMatrix<3,3>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Trace of Matrix
////////////////////////////////////////////////////////////////////////////////

/// Trace of a Matrix
/**
 * Returns the Trace of a Matrix.
 */
/// \{
template <typename T>
inline typename MathMatrix<1,1,T>::value_type
Trace(const MathMatrix<1,1,T>& m);

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
Trace(const MathMatrix<2,2,T>& m);

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
Trace(const MathMatrix<3,3,T>& m);
/// \}

////////////////////////////////////////////////////////////////////////////////
// Scalar operations for Matrices
////////////////////////////////////////////////////////////////////////////////

/// Set each matrix entry to a scalar (componentwise)
template <typename matrix_t>
inline void
MatSet(matrix_t& mInOut, typename matrix_t::value_type s);

/// Set each diagonal of a matrix to a scalar (componentwise)
template <typename matrix_t>
inline void
MatDiagSet(matrix_t& mInOut, typename matrix_t::value_type s);

/// Add a scalar to a matrix (componentwise)
template <typename matrix_t>
inline void
MatAdd(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

/// Subtract a scalar from a matrix (componentwise)
template <typename matrix_t>
inline void
MatSubtract(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

/// Devide a matrix by a scalar (componentwise)
template <typename matrix_t>
inline void
MatDevide(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

/// Multiply a matrix by a scalar (componentwise)
template <typename matrix_t>
inline void
MatMultiply(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s);

///	Fills the matrix with the identity matrix.
template <typename matrix_t>
inline void
MatIdentity(matrix_t& mOut);

///	Fills the matrix with a matrix that rotates around the x-axis in 3 dimensions.
/**	matrix_t has to have at least 3 rows and 3 columns.
 *	Rotation is specified in radiants.*/
template <typename matrix_t>
inline void
MatRotationX(matrix_t& mOut, typename matrix_t::value_type rads);

///	Fills the matrix with a matrix that rotates around the y-axis in 3 dimensions.
/**	matrix_t has to have at least 3 rows and 3 columns.
 *	Rotation is specified in radiants.*/
template <typename matrix_t>
inline void
MatRotationY(matrix_t& mOut, typename matrix_t::value_type rads);

///	Fills the matrix with a matrix that rotates around the y-axis in 2 or 3 dimensions.
/**	matrix_t has to have at least 2 rows and 2 columns.
 *	Rotation is specified in radiants.*/
template <typename matrix_t>
inline void
MatRotationZ(matrix_t& mOut, typename matrix_t::value_type rads);

////////////////////////////////////////////////////////////////////////////////
///	Creates a rotation matrix given yaw, pitch and roll in radiants.
////////////////////////////////////////////////////////////////////////////////
template <typename matrix_t>
inline void
MatRotationYawPitchRoll(matrix_t& mOut,
						typename matrix_t::value_type yaw,
						typename matrix_t::value_type pitch,
						typename matrix_t::value_type roll);

////////////////////////////////////////////////////////////////////////////////
///	Creates a householder matrix given the orthogonal vector to the
///	householder-hypersphere through the origin.
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t, typename vector_t>
inline void
MatHouseholder(matrix_t& mOut, const vector_t& orthoVec);

////////////////////////////////////////////////////////////////////////////////
// Norms for Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline typename matrix_t::value_type
MatFrobeniusNormSq(matrix_t& m);

template <typename matrix_t>
inline typename matrix_t::value_type
MatFrobeniusNorm(matrix_t& m);

template <typename matrix_t>
inline typename matrix_t::value_type
MatOneNorm(matrix_t& m);

template <typename matrix_t>
inline typename matrix_t::value_type
MatInftyNorm(matrix_t& m);

template <typename matrix_t>
inline typename matrix_t::value_type
MatMaxNorm(matrix_t& m);

/// Computes maximum eigenvalue of a (symmetric) matrix
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
MaxAbsEigenvalue(const MathMatrix<M,N,T>& m);

/// Computes minimum eigenvalue of a (symmetric) matrix
template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
MinAbsEigenvalue(const MathMatrix<M,N,T>& m);


// end group math_matrix
/// \}

} //end of namespace

//	include a general, but not very fast implementation of the declared methods above.
#include "math_matrix_functions_common_impl.hpp"

#endif /* __H__UG__COMMON__MATH_MATRIX_FUNCTIONS__ */

/*
 * lgmath_matrix_functions_common_impl.hpp
 *
 *  Created on: 07.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LGMATH__LGMATH_MATRIX_FUNCTIONS_COMMON_IMPL__
#define __H__LGMATH__LGMATH_MATRIX_FUNCTIONS_COMMON_IMPL__

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "math_matrix.h"
#include "common/static_assert.h"

namespace ug
{
///	adds two matrices and stores the result in a third one
// mOut = m1 + m2
template <typename matrix_t>
inline
void
MatAdd(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m1(i, j) + m2(i, j);
		}
}

///	subtracts m2 from m1 and stores the result in a mOut
// mOut = m1 - m2
template <typename matrix_t>
inline
void
MatSubtract(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m1(i, j) - m2(i, j);
		}
}

///	multiply two matrices and stores the result in a third one
// mOut = m1 * m2
template <size_t N, size_t M, size_t L, typename T>
inline
void
MatMultiply(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1, const MathMatrix<L, M, T>& m2)
{
	for(size_t i = 0; i < N; ++i)
		for(size_t j = 0; j < M; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < L; ++k)
			{
				mOut(i,j) += m1(i,k) * m2(k,j);
			}
		}
}

///	multiply two transposed matrices and stores the result in a third one
// mOut = m1^T * m2^T
template <size_t N, size_t M, size_t L, typename T>
inline
void
MatMultiplyTransposed(MathMatrix<N, M, T>& mOut, const MathMatrix<L, N, T>& m1, const MathMatrix<M, L, T>& m2)
{
	for(size_t i = 0; i < N; ++i)
		for(size_t j = 0; j < M; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < L; ++k)
			{
				mOut(i,j) += m1(k,i) * m2(j,k);
			}
		}
}

///	multiply a transposed matrix with itself and stores the result in a second one
// mOut = m^T * m
template <size_t N, size_t M, typename T>
inline
void
MatMultiplyMTM(MathMatrix<N, N, T>& mOut, const MathMatrix<M, N, T>& m)
{
	for(size_t i = 0; i < N; ++i)
	{
		for(size_t j = 0; j < i; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < M; ++k)
			{
				mOut(i,j) += m(k,i) * m(k,j);
			}
			mOut(j,i) = mOut(i,j);
		}
		mOut(i,i) = 0;
		for(size_t k = 0; k < M; ++k)
		{
			mOut(i,i) += m(k,i) * m(k,i);
		}
	}
}

///	multiply a matrix with its transposed and stores the result in a second one
// mOut = m * m^T
template <size_t N, size_t M, typename T>
inline
void
MatMultiplyMMT(MathMatrix<M, M, T>& mOut, const MathMatrix<M, N, T>& m)
{
	for(size_t i = 0; i < M; ++i)
	{
		for(size_t j = 0; j < i; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < N; ++k)
			{
				mOut(i,j) += m(i,k) * m(j,k);
			}
			mOut(j,i) = mOut(i,j);
		}
		mOut(i,i) = 0;
		for(size_t k = 0; k < N; ++k)
		{
			mOut(i,i) += m(i,k) * m(i,k);
		}
	}
}

///	scales a matrix_t and returns the resulting matrix_t
// mOut = scaleFac * m
template <typename matrix_t>
inline
void
MatScale(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i, j) * s;
		}
}

/// transpose a matrix_t
/// transpose a matrix_t
template <typename matrix_t>
inline
void
Transpose(matrix_t& mOut, const matrix_t& m)
{
	assert(&mOut != &m && "ERROR in Transpose(matrix_t& mOut, const matrix_t& m): mOut and m have to be different");

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(j, i);
		}
}

/// transpose a matrix_t, override original matrix_t
template <typename matrix_t>
inline
void
Transpose(matrix_t& m)
{
	assert(m.num_rows()==m.num_cols() && "ERROR in Transpose(matrix_t& m): Square Matrix needed");

	typedef typename matrix_t::size_type size_type;
	matrix_t _temp;
	for(size_type i = 1; i < m.num_rows(); ++i)
		for(size_type j = 0; j < i; ++j)
			_temp(i, j) = m(i, j);

	for(size_type i = 1; i < m.num_rows(); ++i)
		for(size_type j = 0; j < i; ++j)
			m(i, j) = m(j, i);

	for(size_type i = 0; i < m.num_rows()-1; ++i)
		for(size_type j = i+1; j < m.num_cols(); ++j)
			m(i, j) = _temp(j, i);
}

template <typename T>
inline
typename MathMatrix<1,1,T>::value_type
Determinant(const MathMatrix<1,1,T>& m)
{
	return m(0,0);
}

template <typename T>
inline
typename MathMatrix<2,2,T>::value_type
Determinant(const MathMatrix<2,2,T>& m)
{
	return (m(0,0)*m(1,1) - m(1,0)*m(0,1));
}

template <typename T>
inline
typename MathMatrix<3,3,T>::value_type
Determinant(const MathMatrix<3,3,T>& m)
{
	return 	m(0,0)*m(1,1)*m(2,2)
			+ m(0,1)*m(1,2)*m(2,0)
			+ m(0,2)*m(1,0)*m(2,1)
			- m(0,0)*m(1,2)*m(2,1)
			- m(0,1)*m(1,0)*m(2,2)
			- m(0,2)*m(1,1)*m(2,0);
}

template <typename T>
inline
void
Inverse(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m)
{
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(m(0,0) != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	mOut(0,0) =  1./m(0,0);
}

template <typename T>
inline
void
Inverse(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m)
{
	typename MathMatrix<2,2,T>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	typename MathMatrix<2,2,T>::value_type invdet = 1./det;

	mOut(0,0) =  m(1,1)*invdet;
	mOut(1,0) = -m(1,0)*invdet;
	mOut(0,1) = -m(0,1)*invdet;
	mOut(1,1) =  m(0,0)*invdet;
}

template <typename T>
inline
void
Inverse(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m)
{
	typename MathMatrix<3,3,T>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	typename MathMatrix<3,3,T>::value_type invdet = 1./det;

	mOut(0,0) = ( m(1,1)*m(2,2) - m(1,2)*m(2,1)) * invdet;
	mOut(0,1) = (-m(0,1)*m(2,2) + m(0,2)*m(2,1)) * invdet;
	mOut(0,2) = ( m(0,1)*m(1,2) - m(0,2)*m(1,1)) * invdet;
	mOut(1,0) = (-m(1,0)*m(2,2) + m(1,2)*m(2,0)) * invdet;
	mOut(1,1) = ( m(0,0)*m(2,2) - m(0,2)*m(2,0)) * invdet;
	mOut(1,2) = (-m(0,0)*m(1,2) + m(0,2)*m(1,0)) * invdet;
	mOut(2,0) = ( m(1,0)*m(2,1) - m(1,1)*m(2,0)) * invdet;
	mOut(2,1) = (-m(0,0)*m(2,1) + m(0,1)*m(2,0)) * invdet;
	mOut(2,2) = ( m(0,0)*m(1,1) - m(0,1)*m(1,0)) * invdet;
}

template <typename T>
inline
void
Inverse(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m, typename MathMatrix<1,1,T>::value_type& det)
{
	det = m(0,0);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(m(0,0) != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	mOut(0,0) =  1./m(0,0);
}


template <typename T>
inline
void
Inverse(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m, typename MathMatrix<2,2,T>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	typename MathMatrix<2,2,T>::value_type invdet = 1./det;

	mOut(0,0) =  m(1,1)*invdet;
	mOut(1,0) = -m(1,0)*invdet;
	mOut(0,1) = -m(0,1)*invdet;
	mOut(1,1) =  m(0,0)*invdet;
}

template <typename T>
inline
void
Inverse(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m, typename MathMatrix<3,3,T>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	typename MathMatrix<3,3,T>::value_type invdet = 1./det;

	mOut(0,0) = ( m(1,1)*m(2,2) - m(1,2)*m(2,1)) * invdet;
	mOut(0,1) = (-m(0,1)*m(2,2) + m(0,2)*m(2,1)) * invdet;
	mOut(0,2) = ( m(0,1)*m(1,2) - m(0,2)*m(1,1)) * invdet;
	mOut(1,0) = (-m(1,0)*m(2,2) + m(1,2)*m(2,0)) * invdet;
	mOut(1,1) = ( m(0,0)*m(2,2) - m(0,2)*m(2,0)) * invdet;
	mOut(1,2) = (-m(0,0)*m(1,2) + m(0,2)*m(1,0)) * invdet;
	mOut(2,0) = ( m(1,0)*m(2,1) - m(1,1)*m(2,0)) * invdet;
	mOut(2,1) = (-m(0,0)*m(2,1) + m(0,1)*m(2,0)) * invdet;
	mOut(2,2) = ( m(0,0)*m(1,1) - m(0,1)*m(1,0)) * invdet;
}

template <typename T>
inline
void
InverseTransposed(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m)
{
	Inverse(mOut, m);
}

template <typename T>
inline
void
InverseTransposed(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m)
{
	typename MathMatrix<2,2,T>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	typename MathMatrix<2,2,T>::value_type invdet = 1./det;

	mOut(0,0) =  m(1,1)*invdet;
	mOut(1,0) = -m(0,1)*invdet;
	mOut(0,1) = -m(1,0)*invdet;
	mOut(1,1) =  m(0,0)*invdet;
}

template <typename T>
inline
void
InverseTransposed(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m)
{
	typename MathMatrix<3,3,T>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	typename MathMatrix<3,3,T>::value_type invdet = 1./det;

    mOut(0,0) = ( m(1,1)*m(2,2) - m(2,1)*m(1,2)) * invdet;
    mOut(0,1) = (-m(1,0)*m(2,2) + m(2,0)*m(1,2)) * invdet;
    mOut(0,2) = ( m(1,0)*m(2,1) - m(2,0)*m(1,1)) * invdet;
    mOut(1,0) = (-m(0,1)*m(2,2) + m(2,1)*m(0,2)) * invdet;
    mOut(1,1) = ( m(0,0)*m(2,2) - m(2,0)*m(0,2)) * invdet;
    mOut(1,2) = (-m(0,0)*m(2,1) + m(2,0)*m(0,1)) * invdet;
    mOut(2,0) = ( m(0,1)*m(1,2) - m(1,1)*m(0,2)) * invdet;
    mOut(2,1) = (-m(0,0)*m(1,2) + m(1,0)*m(0,2)) * invdet;
    mOut(2,2) = ( m(0,0)*m(1,1) - m(1,0)*m(0,1)) * invdet;
}

template <typename T>
inline
void
InverseTransposed(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m, typename MathMatrix<1,1,T>::value_type& det)
{
	Inverse(mOut, m, det);
}

template <typename T>
inline
void
InverseTransposed(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m, typename MathMatrix<2,2,T>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	typename MathMatrix<2,2,T>::value_type invdet = 1./det;

	mOut(0,0) =  m(1,1)*invdet;
	mOut(1,0) = -m(0,1)*invdet;
	mOut(0,1) = -m(1,0)*invdet;
	mOut(1,1) =  m(0,0)*invdet;
}

template <typename T>
inline
void
InverseTransposed(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m, typename MathMatrix<3,3,T>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	typename MathMatrix<3,3,T>::value_type invdet = 1./det;

    mOut(0,0) = ( m(1,1)*m(2,2) - m(2,1)*m(1,2)) * invdet;
    mOut(0,1) = (-m(1,0)*m(2,2) + m(2,0)*m(1,2)) * invdet;
    mOut(0,2) = ( m(1,0)*m(2,1) - m(2,0)*m(1,1)) * invdet;
    mOut(1,0) = (-m(0,1)*m(2,2) + m(2,1)*m(0,2)) * invdet;
    mOut(1,1) = ( m(0,0)*m(2,2) - m(2,0)*m(0,2)) * invdet;
    mOut(1,2) = (-m(0,0)*m(2,1) + m(2,0)*m(0,1)) * invdet;
    mOut(2,0) = ( m(0,1)*m(1,2) - m(1,1)*m(0,2)) * invdet;
    mOut(2,1) = (-m(0,0)*m(1,2) + m(1,0)*m(0,2)) * invdet;
    mOut(2,2) = ( m(0,0)*m(1,1) - m(1,0)*m(0,1)) * invdet;
}

// Right-Inverse of Matrix
template <size_t N, size_t M, typename T>
inline
void
RightInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m)
{
	UG_STATIC_ASSERT(M <= N, pseudo_inverse_does_not_exist);

	// H = m*mT (H is symmetric)
	// TODO: Since H is symmetric, we could store only lower or upper elements
	MathMatrix<M,M,T> H, HInv;
	MatMultiplyMMT(H, m);
	// Invert H
	Inverse(HInv, H);

	MatMultiplyTransposed(mOut, m, HInv);
}

// Left-Inverse of Matrix
template <size_t N, size_t M, typename T>
inline
void
LeftInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m)
{
	UG_STATIC_ASSERT(N <= M, pseudo_inverse_does_not_exist);

	// H = mT*m (H is symmetric)
	// TODO: Since H is symmetric, we could store only lower or upper elements
	MathMatrix<M,M,T> H, HInv;
	MatMultiplyMTM(H, m);

	// Invert H
	Inverse(HInv, H);

	MatMultiplyTransposed(mOut, HInv, m);
}

template <>
inline void
RightInverse(MathMatrix<1,1>& mOut, MathMatrix<1,1>& m){Inverse(mOut, m);}

template <>
inline void
RightInverse(MathMatrix<2,2>& mOut, MathMatrix<2,2>& m){Inverse(mOut, m);}

template <>
inline void
RightInverse(MathMatrix<3,3>& mOut, MathMatrix<3,3>& m){Inverse(mOut, m);}

template <>
inline void
LeftInverse(MathMatrix<1,1>& mOut, MathMatrix<1,1>& m){Inverse(mOut, m);}

template <>
inline void
LeftInverse(MathMatrix<2,2>& mOut, MathMatrix<2,2>& m){Inverse(mOut, m);}

template <>
inline void
LeftInverse(MathMatrix<3,3>& mOut, MathMatrix<3,3>& m){Inverse(mOut, m);}

/// Set each matrix entry to a scalar (componentwise)
template <typename matrix_t>
inline
void
MatSet(matrix_t& mInOut, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mInOut.num_rows(); ++i)
		for(size_type j = 0; j < mInOut.num_cols(); ++j)
		{
			mInOut(i, j) = s;
		}
}

/// Set each diagonal of a matrix to a scalar (componentwise)
template <typename matrix_t>
inline
void
MatDiagSet(matrix_t& mInOut, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mInOut.num_rows(); ++i)
		for(size_type j = 0; j < mInOut.num_cols(); ++j)
		{
			mInOut(i, i) = s;
		}
}

/// Add a scalar to a vector (componentwise)
template <typename matrix_t>
inline
void
MatAdd(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) + s;
		}
}

/// Subtract a scalar from a vector (componentwise)
template <typename matrix_t>
inline
void
MatSubtract(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) - s;
		}
}

/// Devide a vector by a scalar (componentwise)
template <typename matrix_t>
inline
void
MatDevide(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) /s;
		}
}

/// Multiply a vector by a scalar (componentwise)
template <typename matrix_t>
inline
void
MatMultiply(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) * s;
		}
}

template <typename matrix_t>
inline
void
MatIdentity(matrix_t& mOut)
{
	MatSet(mOut, 0);
	MatDiagSet(mOut, 1);
}

template <typename matrix_t>
inline
void
MatRotationX(matrix_t& mOut, typename matrix_t::value_type rads)
{
	//UG_STATIC_ASSERT(matrix_t::RowSize > 2, AT_LEAST_3_ROWS_REQUIRED);
	//UG_STATIC_ASSERT(matrix_t::ColSize > 2, AT_LEAST_3_COLS_REQUIRED);
	
	MatIdentity(mOut);
	typename matrix_t::value_type s = sin(rads);
	typename matrix_t::value_type c = cos(rads);
	
	mOut(1, 1) = c;		mOut(1, 2) = -s;
	mOut(2, 1) = s;		mOut(2, 2) = c;
}

template <typename matrix_t>
inline
void
MatRotationY(matrix_t& mOut, typename matrix_t::value_type rads)
{
	//UG_STATIC_ASSERT(matrix_t::RowSize > 2, AT_LEAST_3_ROWS_REQUIRED);
	//UG_STATIC_ASSERT(matrix_t::ColSize > 2, AT_LEAST_3_COLS_REQUIRED);
	
	MatIdentity(mOut);
	typename matrix_t::value_type s = sin(rads);
	typename matrix_t::value_type c = cos(rads);
	
	mOut(0, 0) = c;			mOut(0, 2) = s;
	mOut(2, 0) = -s;		mOut(2, 2) = c;
}

template <typename matrix_t>
inline
void
MatRotationZ(matrix_t& mOut, typename matrix_t::value_type rads)
{
	MatIdentity(mOut);
	typename matrix_t::value_type s = sin(rads);
	typename matrix_t::value_type c = cos(rads);
	
	mOut(0, 0) = c;		mOut(0, 1) = -s;
	mOut(1, 0) = s;		mOut(1, 1) = c;
}

template <typename matrix_t>
inline
void
MatRotationYawPitchRoll(matrix_t& mOut,
						typename matrix_t::value_type yaw,
						typename matrix_t::value_type pitch,
						typename matrix_t::value_type roll)
{
	//UG_STATIC_ASSERT(matrix_t::RowSize > 2, AT_LEAST_3_ROWS_REQUIRED);
	//UG_STATIC_ASSERT(matrix_t::ColSize > 2, AT_LEAST_3_COLS_REQUIRED);

	matrix_t tMat1, tMat2, tMat3;
	MatRotationX(tMat1, yaw);
	MatRotationY(tMat2, pitch);
	MatMultiply(tMat3, tMat1, tMat2);
	MatRotationZ(tMat1, roll);
	MatMultiply(mOut, tMat3, tMat1);
}

template <typename matrix_t>
inline
typename matrix_t::value_type
MatFrobeniusNormSq(matrix_t& m)
{
	typename matrix_t::value_type norm = 0;
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < m.num_rows(); ++i)
		for(size_type j = 0; j < m.num_cols(); ++j)
		{
			norm += m(i,j)*m(i,j);
		}

	return norm;
}

template <typename matrix_t>
inline
typename matrix_t::value_type
MatFrobeniusNorm(matrix_t& m)
{
	return static_cast<typename matrix_t::value_type>(sqrt(MatFrobeniusNormSq(m)));
}

template <typename matrix_t>
inline
typename matrix_t::value_type
MatOneNorm(matrix_t& m)
{
	typename matrix_t::value_type sum, max = 0;
	typedef typename matrix_t::size_type size_type;
	for(size_type j = 0; j < m.num_cols(); ++j)
	{
		sum = 0;
		for(size_type i = 0; i < m.num_rows(); ++i)
		{
			sum += fabs(m(i,j));
		}
	max = (sum > max) ? sum : max;
	}
	return max;
}

template <typename matrix_t>
inline
typename matrix_t::value_type
MatInftyNorm(matrix_t& m)
{
	typename matrix_t::value_type sum, max = 0;
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < m.num_rows(); ++i)
	{
		sum = 0;
		for(size_type j = 0; j < m.num_cols(); ++j)
		{
			sum += fabs(m(i,j));
		}
	max = (sum > max) ? sum : max;
	}
	return max;
}

template <typename matrix_t>
inline
typename matrix_t::value_type
MatMaxNorm(matrix_t& m)
{
	typename matrix_t::value_type max = 0;
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < m.num_rows(); ++i)
		for(size_type j = 0; j < m.num_cols(); ++j)
		{
			max = (m(i,j) > max) ? m(i,j) : max;
		}

	return max;
}

} // end of namespace: lgmath

#endif /* __H__LGMATH__LGMATH_MATRIX_FUNCTIONS_COMMON_IMPL__ */

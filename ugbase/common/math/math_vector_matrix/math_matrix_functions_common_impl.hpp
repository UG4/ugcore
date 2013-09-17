/*
 * math_matrix_functions_common_impl.hpp
 *
 *  Created on: 07.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__COMMON__MATH_MATRIX_FUNCTIONS_COMMON_IMPL__
#define __H__UG__COMMON__MATH_MATRIX_FUNCTIONS_COMMON_IMPL__

#include <cmath>
#include <iostream>
#include <iomanip>
#include <cassert>
#include "math_matrix.h"
#include "common/assert.h"
#include "common/static_assert.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
// Addition of Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline void
MatAdd(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i,j) = m1(i,j) + m2(i,j);
		}
}

////////////////////////////////////////////////////////////////////////////////
// Subtraction of Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline void
MatSubtract(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i,j) = m1(i,j) - m2(i,j);
		}
}


////////////////////////////////////////////////////////////////////////////////
// Multiplication of Matrices
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiply(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1,
            const MathMatrix<L, M, T>& m2)
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

template <size_t N, size_t M, size_t L, size_t P, typename T>
inline void
MatMultiply(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1,
            const MathMatrix<L, P, T>& m2, const MathMatrix<P, M, T>& m3)
{
	MathMatrix<L, M, T> help;

	for(size_t i = 0; i < N; ++i)
		for(size_t j = 0; j < M; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < L; ++k)
			{
				help(k,j) = 0;
				for(size_t l = 0; l < P; ++l)
				{
					help(k,j) += m2(k,l) * m3(l,j);
				}

				mOut(i,j) += m1(i,k) * help(k,j);
			}
		}
}

template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiplyTransposed(MathMatrix<N, M, T>& mOut, const MathMatrix<L, N, T>& m1,
                      const MathMatrix<M, L, T>& m2)
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

template <size_t N, size_t M, typename T>
inline void
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

template <size_t N, size_t M, typename T>
inline void
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

////////////////////////////////////////////////////////////////////////////////
// "Contraction" for Matrices (note: contraction is only known regarding tensors!)
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline typename matrix_t::value_type
MatContraction(const matrix_t& m1, const matrix_t& m2)
{
	typename matrix_t::value_type norm = 0;
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < m1.num_rows(); ++i)
		for(size_type j = 0; j < m1.num_cols(); ++j)
		{
			norm += m1(i,j)*m2(i,j);
		}

	return norm;
}

////////////////////////////////////////////////////////////////////////////////
// Scaling of Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline void
MatScale(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i, j) * s;
		}
}

template <typename matrix_t>
inline void
MatScaleAppend(matrix_t& mOut, typename matrix_t::value_type s, const matrix_t& m)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) += m(i, j) * s;
		}
}

////////////////////////////////////////////////////////////////////////////////
// Transposed of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline void
Transpose(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m)
{
	typedef typename MathMatrix<N,M,T>::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(j, i);
		}
}

template <typename matrix_t>
inline void
Transpose(matrix_t& m)
{
	UG_ASSERT(m.num_rows()==m.num_cols(), "Transpose: Square Matrix needed");

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

////////////////////////////////////////////////////////////////////////////////
// Determinant of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, typename T>
inline typename MathMatrix<N,N,T>::value_type
Determinant(const MathMatrix<N,N,T>& m)
{
	UG_THROW("Determinant for matrix of size "<<N<<"x"<<N<<" not implemented.");
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
Determinant(const MathMatrix<1,1,T>& m)
{
	return m(0,0);
}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
Determinant(const MathMatrix<2,2,T>& m)
{
	return (m(0,0)*m(1,1) - m(1,0)*m(0,1));
}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
Determinant(const MathMatrix<3,3,T>& m)
{
	return 	m(0,0)*m(1,1)*m(2,2)
			+ m(0,1)*m(1,2)*m(2,0)
			+ m(0,2)*m(1,0)*m(2,1)
			- m(0,0)*m(1,2)*m(2,1)
			- m(0,1)*m(1,0)*m(2,2)
			- m(0,2)*m(1,1)*m(2,0);
}

////////////////////////////////////////////////////////////////////////////////
// Gram Determinant of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
GramDeterminant(const MathMatrix<N,M,T>& m)
{
	if(N <= M)
	{
		MathMatrix<N,N,T> mmT;
		MatMultiplyMMT(mmT, m);
		return Determinant(mmT);
	}
	else
	{
		MathMatrix<M,M,T> mTm;
		MatMultiplyMTM(mTm, m);
		return Determinant(mTm);
	}
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
GramDeterminant(const MathMatrix<1,1,T>& m)
{
	return pow(Determinant(m), 2);
}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
GramDeterminant(const MathMatrix<2,2,T>& m)
{
	return pow(Determinant(m), 2);
}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
GramDeterminant(const MathMatrix<3,3,T>& m)
{
	return pow(Determinant(m), 2);
}

////////////////////////////////////////////////////////////////////////////////
// Sqrt Gram Determinant of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
SqrtGramDeterminant(const MathMatrix<N,M,T>& m)
{
	return sqrt(GramDeterminant(m));
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
SqrtGramDeterminant(const MathMatrix<1,1,T>& m)
{
	return fabs(Determinant(m));
}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
SqrtGramDeterminant(const MathMatrix<2,2,T>& m)
{
	return fabs(Determinant(m));
}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
SqrtGramDeterminant(const MathMatrix<3,3,T>& m)
{
	return fabs(Determinant(m));
}

////////////////////////////////////////////////////////////////////////////////
// Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
Inverse(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m)
{
	UG_THROW("Inverse for matrix of size "<<N<<"x"<<M<<" not implemented.");
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
Inverse(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m)
{
	const typename MathMatrix<1,1,T>::value_type det = m(0,0);
	UG_ASSERT(&mOut != &m, "Inverse: mOut and m have to be different");
	UG_ASSERT(det != 0, "Inverse: determinate is zero, can not Invert Matrix");
	mOut(0,0) =  1./m(0,0);
	return det;
}


template <typename T>
inline typename MathMatrix<2,2,T>::value_type
Inverse(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m)
{
	const typename MathMatrix<2,2,T>::value_type det = Determinant(m);
	UG_ASSERT(&mOut != &m, "Inverse: mOut and m have to be different");
	UG_ASSERT(det != 0, "Inverse: determinate is zero, can not Invert Matrix");
	const typename MathMatrix<2,2,T>::value_type invdet = 1./det;

	mOut(0,0) =  m(1,1)*invdet;
	mOut(1,0) = -m(1,0)*invdet;
	mOut(0,1) = -m(0,1)*invdet;
	mOut(1,1) =  m(0,0)*invdet;

	return det;
}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
Inverse(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m)
{
	const typename MathMatrix<3,3,T>::value_type det = Determinant(m);
	UG_ASSERT(&mOut != &m, "Inverse: mOut and m have to be different");
	UG_ASSERT(det != 0, "Inverse: determinate is zero, can not Invert Matrix");
	const typename MathMatrix<3,3,T>::value_type invdet = 1./det;

	mOut(0,0) = ( m(1,1)*m(2,2) - m(1,2)*m(2,1)) * invdet;
	mOut(0,1) = (-m(0,1)*m(2,2) + m(0,2)*m(2,1)) * invdet;
	mOut(0,2) = ( m(0,1)*m(1,2) - m(0,2)*m(1,1)) * invdet;
	mOut(1,0) = (-m(1,0)*m(2,2) + m(1,2)*m(2,0)) * invdet;
	mOut(1,1) = ( m(0,0)*m(2,2) - m(0,2)*m(2,0)) * invdet;
	mOut(1,2) = (-m(0,0)*m(1,2) + m(0,2)*m(1,0)) * invdet;
	mOut(2,0) = ( m(1,0)*m(2,1) - m(1,1)*m(2,0)) * invdet;
	mOut(2,1) = (-m(0,0)*m(2,1) + m(0,1)*m(2,0)) * invdet;
	mOut(2,2) = ( m(0,0)*m(1,1) - m(0,1)*m(1,0)) * invdet;

	return det;
}

////////////////////////////////////////////////////////////////////////////////
// Inverse Transposed of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
InverseTransposed(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m)
{
	UG_THROW("InverseTransposed for matrix of size "<<M<<"x"<<N<<" not implemented.");
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
InverseTransposed(MathMatrix<1,1,T>& mOut, const MathMatrix<1,1,T>& m)
{
	return Inverse(mOut, m);
}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
InverseTransposed(MathMatrix<2,2,T>& mOut, const MathMatrix<2,2,T>& m)
{
	const typename MathMatrix<2,2,T>::value_type det = Determinant(m);
	UG_ASSERT(&mOut != &m, "InverseTransposed: mOut and m have to be different");
	UG_ASSERT(det != 0, "InverseTransposed: determinate is zero, can not Invert Matrix");
	const typename MathMatrix<2,2,T>::value_type invdet = 1./det;

	mOut(0,0) =  m(1,1)*invdet;
	mOut(1,0) = -m(0,1)*invdet;
	mOut(0,1) = -m(1,0)*invdet;
	mOut(1,1) =  m(0,0)*invdet;

	return det;
}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
InverseTransposed(MathMatrix<3,3,T>& mOut, const MathMatrix<3,3,T>& m)
{
	const typename MathMatrix<3,3,T>::value_type det = Determinant(m);
	UG_ASSERT(&mOut != &m, "InverseTransposed: mOut and m have to be different");
	UG_ASSERT(det != 0, "InverseTransposed: determinate is zero, can not Invert Matrix");
	const typename MathMatrix<3,3,T>::value_type invdet = 1./det;

    mOut(0,0) = ( m(1,1)*m(2,2) - m(2,1)*m(1,2)) * invdet;
    mOut(0,1) = (-m(1,0)*m(2,2) + m(2,0)*m(1,2)) * invdet;
    mOut(0,2) = ( m(1,0)*m(2,1) - m(2,0)*m(1,1)) * invdet;
    mOut(1,0) = (-m(0,1)*m(2,2) + m(2,1)*m(0,2)) * invdet;
    mOut(1,1) = ( m(0,0)*m(2,2) - m(2,0)*m(0,2)) * invdet;
    mOut(1,2) = (-m(0,0)*m(2,1) + m(2,0)*m(0,1)) * invdet;
    mOut(2,0) = ( m(0,1)*m(1,2) - m(1,1)*m(0,2)) * invdet;
    mOut(2,1) = (-m(0,0)*m(1,2) + m(1,0)*m(0,2)) * invdet;
    mOut(2,2) = ( m(0,0)*m(1,1) - m(1,0)*m(0,1)) * invdet;

    return det;
}

////////////////////////////////////////////////////////////////////////////////
// Right-Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
RightInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m)
{
	UG_STATIC_ASSERT(M <= N, pseudo_inverse_does_not_exist);

	// H = m*mT (H is symmetric)
	// TODO: Since H is symmetric, we could store only lower or upper elements
	MathMatrix<M,M,T> H, HInv;
	MatMultiplyMMT(H, m);
	// Invert H
	const number det = Inverse(HInv, H);

	MatMultiplyTransposed(mOut, m, HInv);

	return sqrt(det);
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
RightInverse(MathMatrix<1,1>& mOut, MathMatrix<1,1>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
RightInverse(MathMatrix<2,2>& mOut, MathMatrix<2,2>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
RightInverse(MathMatrix<3,3>& mOut, MathMatrix<3,3>& m){return fabs(Inverse(mOut, m));}

////////////////////////////////////////////////////////////////////////////////
// Left-Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
LeftInverse(MathMatrix<N,M,T>& mOut, MathMatrix<M,N,T>& m)
{
	UG_STATIC_ASSERT(N <= M, pseudo_inverse_does_not_exist);

	// H = mT*m (H is symmetric)
	// TODO: Since H is symmetric, we could store only lower or upper elements
	MathMatrix<N,N,T> H, HInv;
	MatMultiplyMTM(H, m);

	// Invert H
	const number det = Inverse(HInv, H);

	MatMultiplyTransposed(mOut, HInv, m);

	return sqrt(det);
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
LeftInverse(MathMatrix<1,1>& mOut, MathMatrix<1,1>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
LeftInverse(MathMatrix<2,2>& mOut, MathMatrix<2,2>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
LeftInverse(MathMatrix<3,3>& mOut, MathMatrix<3,3>& m){return fabs(Inverse(mOut, m));}

////////////////////////////////////////////////////////////////////////////////
// Trace of Matrix
////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
Trace(const MathMatrix<1,1,T>& m)
{
	return m(0,0);
}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
Trace(const MathMatrix<2,2,T>& m)
{
	return (m(0,0)+m(1,1));
}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
Trace(const MathMatrix<3,3,T>& m)
{
	return 	(m(0,0)+m(1,1)+m(2,2));
}

////////////////////////////////////////////////////////////////////////////////
// Scalar operations for Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline void
MatSet(matrix_t& mInOut, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mInOut.num_rows(); ++i)
		for(size_type j = 0; j < mInOut.num_cols(); ++j)
		{
			mInOut(i, j) = s;
		}
}

template <typename matrix_t>
inline void
MatDiagSet(matrix_t& mInOut, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mInOut.num_rows(); ++i)
		for(size_type j = 0; j < mInOut.num_cols(); ++j)
		{
			mInOut(i, i) = s;
		}
}

template <typename matrix_t>
inline void
MatAdd(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) + s;
		}
}

template <typename matrix_t>
inline void
MatSubtract(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) - s;
		}
}

template <typename matrix_t>
inline void
MatDevide(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.num_rows(); ++i)
		for(size_type j = 0; j < mOut.num_cols(); ++j)
		{
			mOut(i, j) = m(i,j) /s;
		}
}

template <typename matrix_t>
inline void
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
inline void
MatIdentity(matrix_t& mOut)
{
	MatSet(mOut, 0);
	MatDiagSet(mOut, 1);
}

template <typename matrix_t>
inline void
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
inline void
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
inline void
MatRotationZ(matrix_t& mOut, typename matrix_t::value_type rads)
{
	MatIdentity(mOut);
	typename matrix_t::value_type s = sin(rads);
	typename matrix_t::value_type c = cos(rads);
	
	mOut(0, 0) = c;		mOut(0, 1) = -s;
	mOut(1, 0) = s;		mOut(1, 1) = c;
}

////////////////////////////////////////////////////////////////////////////////
///	Creates a rotation matrix given yaw, pitch and roll in radiants.
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline void
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

////////////////////////////////////////////////////////////////////////////////
// Norms for Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline typename matrix_t::value_type
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
inline typename matrix_t::value_type
MatFrobeniusNorm(matrix_t& m)
{
	return static_cast<typename matrix_t::value_type>(sqrt(MatFrobeniusNormSq(m)));
}

template <typename matrix_t>
inline typename matrix_t::value_type
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
inline typename matrix_t::value_type
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
inline typename matrix_t::value_type
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



template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
MaxAbsEigenvalue(const MathMatrix<M,N,T>& m)
{
	UG_THROW("MaxAbsEigenvalue for matrix of size "<<N<<"x"<<M<<" not implemented.");
}

template <typename T>
inline typename MathMatrix<1,1,T>::value_type
MaxAbsEigenvalue(const MathMatrix<1,1,T>& m)
{
	const typename MathMatrix<1,1,T>::value_type val = m(0,0);
	return fabs(val);
}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
MaxAbsEigenvalue(const MathMatrix<2,2,T>& m)
{
	typename MathMatrix<2,2,T>::value_type minus_p_half, val;
	minus_p_half = m(0,0)+m(1,1);
	val = minus_p_half*minus_p_half - (m(0,0)*m(1,1) - m(0,1)*m(1,0));
	UG_ASSERT(val >= 0.0, "MaxAbsEigenvalues: Complex Eigenvalues???");

	if (minus_p_half >=0.0) {return (minus_p_half + sqrt(val));}
	else {return fabs(minus_p_half-sqrt(val));}
}


template <typename matrix_t>
inline typename matrix_t::value_type
MinAbsEigenvalue(matrix_t& m)
{
	matrix_t inv;
	Inverse(inv, m);
	typename matrix_t::value_type min=1.0/MaxAbsEigenvalue(inv);
	return min;
}

} // end of namespace

#endif /* __H__UG__COMMON__MATH_MATRIX_FUNCTIONS_COMMON_IMPL__ */

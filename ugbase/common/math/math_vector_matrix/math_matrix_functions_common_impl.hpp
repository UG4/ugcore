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

namespace ug
{
///	adds two matrices and stores the result in a third one
// mOut = m1 + m2
template <int N, int M>
inline
void
MatAdd(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m1, const MathMatrix<N,M>& m2)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m1.entry(i, j) + m2.entry(i, j);
}

///	subtracts m2 from m1 and stores the result in a mOut
// mOut = m1 - m2
template <int N, int M>
inline
void
MatSubtract(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m1, const MathMatrix<N,M>& m2)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m1.entry(i, j) - m2.entry(i, j);
}

///	scales a MathMatrix<N,M> and returns the resulting MathMatrix<N,M>
// mOut = scaleFac * m
template <int N, int M>
inline
void
MatScale(MathMatrix<N,M>& mOut, typename MathMatrix<N,M>::value_type s, const MathMatrix<N,M>& m)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m.entry(i, j) * s;
}

/// transpose a MathMatrix<N,M>
template <int N, int M>
inline
void
Transpose(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m)
{
	assert(&mOut != &m && "ERROR in Transpose(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m): mOut and m have to be different");
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m.entry(j, i);
}

/// transpose a MathMatrix<N,M>, override original MathMatrix<N,M>
template <int N, int M>
inline
void
Transpose(MathMatrix<N,M>& m)
{
	MathMatrix<N,M> _temp;
	for(uint i = 1; i < N; ++i)
		for(uint j = 0; j < i; ++j)
			_temp.entry(i, j) = m.entry(i, j);

	for(uint i = 1; i < N; ++i)
		for(uint j = 0; j < i; ++j)
			m.entry(i, j) = m.entry(j, i);

	for(uint i = 0; i < N-1; ++i)
		for(uint j = i+1; j < M; ++j)
			m.entry(i, j) = _temp.entry(j, i);
}

inline
MathMatrix<2,2>::value_type
Determinant(const MathMatrix<2,2>& m)
{
	return (m.entry(0,0)*m.entry(1,1) - m.entry(1,0)*m.entry(0,1));
}

inline
MathMatrix<3,3>::value_type
Determinant(const MathMatrix<3,3>& m)
{
	return 	m.entry(0,0)*m.entry(1,1)*m.entry(2,2)
			+ m.entry(0,1)*m.entry(1,2)*m.entry(2,0)
			+ m.entry(0,2)*m.entry(1,0)*m.entry(2,1)
			- m.entry(0,0)*m.entry(1,2)*m.entry(2,1)
			- m.entry(0,1)*m.entry(1,0)*m.entry(2,2)
			- m.entry(0,2)*m.entry(1,1)*m.entry(2,0);
}

inline
void
Inverse(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m)
{
	MathMatrix<2,2>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	MathMatrix<2,2>::value_type invdet = 1./det;

	mOut.entry(0,0) =  m.entry(1,1)*invdet;
	mOut.entry(1,0) = -m.entry(1,0)*invdet;
	mOut.entry(0,1) = -m.entry(0,1)*invdet;
	mOut.entry(1,1) =  m.entry(0,0)*invdet;
}

inline
void
Inverse(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m)
{
	MathMatrix<3,3>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	MathMatrix<3,3>::value_type invdet = 1./det;

	mOut.entry(0,0) = ( m.entry(1,1)*m.entry(2,2) - m.entry(1,2)*m.entry(2,1)) * invdet;
	mOut.entry(0,1) = (-m.entry(0,1)*m.entry(2,2) + m.entry(0,2)*m.entry(2,1)) * invdet;
	mOut.entry(0,2) = ( m.entry(0,1)*m.entry(1,2) - m.entry(0,2)*m.entry(1,1)) * invdet;
	mOut.entry(1,0) = (-m.entry(1,0)*m.entry(2,2) + m.entry(1,2)*m.entry(2,0)) * invdet;
	mOut.entry(1,1) = ( m.entry(0,0)*m.entry(2,2) - m.entry(0,2)*m.entry(2,0)) * invdet;
	mOut.entry(1,2) = (-m.entry(0,0)*m.entry(1,2) + m.entry(0,2)*m.entry(1,0)) * invdet;
	mOut.entry(2,0) = ( m.entry(1,0)*m.entry(2,1) - m.entry(1,1)*m.entry(2,0)) * invdet;
	mOut.entry(2,1) = (-m.entry(0,0)*m.entry(2,1) + m.entry(0,1)*m.entry(2,0)) * invdet;
	mOut.entry(2,2) = ( m.entry(0,0)*m.entry(1,1) - m.entry(0,1)*m.entry(1,0)) * invdet;
}

inline
void
Inverse(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m, MathMatrix<2,2>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	MathMatrix<2,2>::value_type invdet = 1./det;

	mOut.entry(0,0) =  m.entry(1,1)*invdet;
	mOut.entry(1,0) = -m.entry(1,0)*invdet;
	mOut.entry(0,1) = -m.entry(0,1)*invdet;
	mOut.entry(1,1) =  m.entry(0,0)*invdet;
}

void
Inverse(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m, MathMatrix<3,3>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in Inverse: determinate is zero, can not Invert Matrix");
	MathMatrix<3,3>::value_type invdet = 1./det;

	mOut.entry(0,0) = ( m.entry(1,1)*m.entry(2,2) - m.entry(1,2)*m.entry(2,1)) * invdet;
	mOut.entry(0,1) = (-m.entry(0,1)*m.entry(2,2) + m.entry(0,2)*m.entry(2,1)) * invdet;
	mOut.entry(0,2) = ( m.entry(0,1)*m.entry(1,2) - m.entry(0,2)*m.entry(1,1)) * invdet;
	mOut.entry(1,0) = (-m.entry(1,0)*m.entry(2,2) + m.entry(1,2)*m.entry(2,0)) * invdet;
	mOut.entry(1,1) = ( m.entry(0,0)*m.entry(2,2) - m.entry(0,2)*m.entry(2,0)) * invdet;
	mOut.entry(1,2) = (-m.entry(0,0)*m.entry(1,2) + m.entry(0,2)*m.entry(1,0)) * invdet;
	mOut.entry(2,0) = ( m.entry(1,0)*m.entry(2,1) - m.entry(1,1)*m.entry(2,0)) * invdet;
	mOut.entry(2,1) = (-m.entry(0,0)*m.entry(2,1) + m.entry(0,1)*m.entry(2,0)) * invdet;
	mOut.entry(2,2) = ( m.entry(0,0)*m.entry(1,1) - m.entry(0,1)*m.entry(1,0)) * invdet;
}

inline
void
InverseTransposed(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m)
{
	MathMatrix<2,2>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	MathMatrix<2,2>::value_type invdet = 1./det;

	mOut.entry(0,0) =  m.entry(1,1)*invdet;
	mOut.entry(1,0) = -m.entry(0,1)*invdet;
	mOut.entry(0,1) = -m.entry(1,0)*invdet;
	mOut.entry(1,1) =  m.entry(0,0)*invdet;
}

inline
void
InverseTransposed(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m)
{
	MathMatrix<3,3>::value_type det;
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	MathMatrix<3,3>::value_type invdet = 1./det;

    mOut.entry(0,0) = ( m.entry(1,1)*m.entry(2,2) - m.entry(2,1)*m.entry(1,2)) * invdet;
    mOut.entry(0,1) = (-m.entry(1,0)*m.entry(2,2) + m.entry(2,0)*m.entry(1,2)) * invdet;
    mOut.entry(0,2) = ( m.entry(1,0)*m.entry(2,1) - m.entry(2,0)*m.entry(1,1)) * invdet;
    mOut.entry(1,0) = (-m.entry(0,1)*m.entry(2,2) + m.entry(2,1)*m.entry(0,2)) * invdet;
    mOut.entry(1,1) = ( m.entry(0,0)*m.entry(2,2) - m.entry(2,0)*m.entry(0,2)) * invdet;
    mOut.entry(1,2) = (-m.entry(0,0)*m.entry(2,1) + m.entry(2,0)*m.entry(0,1)) * invdet;
    mOut.entry(2,0) = ( m.entry(0,1)*m.entry(1,2) - m.entry(1,1)*m.entry(0,2)) * invdet;
    mOut.entry(2,1) = (-m.entry(0,0)*m.entry(1,2) + m.entry(1,0)*m.entry(0,2)) * invdet;
    mOut.entry(2,2) = ( m.entry(0,0)*m.entry(1,1) - m.entry(1,0)*m.entry(0,1)) * invdet;
}

inline
void
InverseTransposed(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m, MathMatrix<2,2>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	MathMatrix<2,2>::value_type invdet = 1./det;

	mOut.entry(0,0) =  m.entry(1,1)*invdet;
	mOut.entry(1,0) = -m.entry(0,1)*invdet;
	mOut.entry(0,1) = -m.entry(1,0)*invdet;
	mOut.entry(1,1) =  m.entry(0,0)*invdet;
}

inline
void
InverseTransposed(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m, MathMatrix<3,3>::value_type& det)
{
	det = Determinant(m);
	assert(&mOut != &m && "ERROR: mOut and m have to be different");
	assert(det != 0 && "ERROR in InverseTransposed: determinate is zero, can not Invert Matrix");
	MathMatrix<3,3>::value_type invdet = 1./det;

    mOut.entry(0,0) = ( m.entry(1,1)*m.entry(2,2) - m.entry(2,1)*m.entry(1,2)) * invdet;
    mOut.entry(0,1) = (-m.entry(1,0)*m.entry(2,2) + m.entry(2,0)*m.entry(1,2)) * invdet;
    mOut.entry(0,2) = ( m.entry(1,0)*m.entry(2,1) - m.entry(2,0)*m.entry(1,1)) * invdet;
    mOut.entry(1,0) = (-m.entry(0,1)*m.entry(2,2) + m.entry(2,1)*m.entry(0,2)) * invdet;
    mOut.entry(1,1) = ( m.entry(0,0)*m.entry(2,2) - m.entry(2,0)*m.entry(0,2)) * invdet;
    mOut.entry(1,2) = (-m.entry(0,0)*m.entry(2,1) + m.entry(2,0)*m.entry(0,1)) * invdet;
    mOut.entry(2,0) = ( m.entry(0,1)*m.entry(1,2) - m.entry(1,1)*m.entry(0,2)) * invdet;
    mOut.entry(2,1) = (-m.entry(0,0)*m.entry(1,2) + m.entry(1,0)*m.entry(0,2)) * invdet;
    mOut.entry(2,2) = ( m.entry(0,0)*m.entry(1,1) - m.entry(1,0)*m.entry(0,1)) * invdet;
}

/// Set each matrix entry to a scalar (componentwise)
template <int N, int M>
inline
void
MatSet(MathMatrix<N,M>& mInOut, typename MathMatrix<N,M>::value_type s)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mInOut.entry(i, j) = s;
}

/// Set each diagonal of a matrix to a scalar (componentwise)
template <int N, int M>
inline
void
MatDiagSet(MathMatrix<N,M>& mInOut, typename MathMatrix<N,M>::value_type s)
{
	for(uint i = 0; i < N; ++i)
		mInOut.entry(i, i) = s;
}

/// Add a scalar to a vector (componentwise)
template <int N, int M>
inline
void
MatAdd(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m.entry(i,j) + s;
}

/// Subtract a scalar from a vector (componentwise)
template <int N, int M>
inline
void
MatSubtract(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m.entry(i,j) - s;
}

/// Devide a vector by a scalar (componentwise)
template <int N, int M>
inline
void
MatDevide(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m.entry(i,j) /s;
}

/// Multiply a vector by a scalar (componentwise)
template <int N, int M>
inline
void
MatMultiply(MathMatrix<N,M>& mOut, const MathMatrix<N,M>& m, typename MathMatrix<N,M>::value_type s)
{
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			mOut.entry(i, j) = m.entry(i,j) * s;
}

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatFrobeniusNormSq(MathMatrix<N,M>& m)
{
	typename MathMatrix<N,M>::value_type norm = 0;
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			norm += m.entry(i,j)*m.entry(i,j);

	return norm;
}

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatFrobeniusNorm(MathMatrix<N,M>& m)
{
	return static_cast<typename MathMatrix<N,M>::value_type>(sqrt(MatFrobeniusNormSq(m)));
}

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatOneNorm(MathMatrix<N,M>& m)
{
	typename MathMatrix<N,M>::value_type sum, max = 0;
	for(uint j = 0; j < N; ++j)
	{
		sum = 0;
		for(uint i = 0; i < M; ++i)
		{
			sum += fabs(m.entry(i,j));
		}
	max = (sum > max) ? sum : max;
	}
	return max;
}

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatInftyNorm(MathMatrix<N,M>& m)
{
	typename MathMatrix<N,M>::value_type sum, max = 0;
	for(uint i = 0; i < N; ++i)
	{
		sum = 0;
		for(uint j = 0; j < M; ++j)
		{
			sum += fabs(m.entry(i,j));
		}
	max = (sum > max) ? sum : max;
	}
	return max;
}

template <int N, int M>
inline
typename MathMatrix<N,M>::value_type
MatMaxNorm(MathMatrix<N,M>& m)
{
	typename MathMatrix<N,M>::value_type max = 0;
	for(uint i = 0; i < N; ++i)
		for(uint j = 0; j < M; ++j)
			max = (m.entry(i,j) > max) ? m.entry(i,j) : max;

	return max;
}

} // end of namespace: lgmath

#endif /* __H__LGMATH__LGMATH_MathMatrix<N,M>_FUNCTIONS_COMMON_IMPL__ */

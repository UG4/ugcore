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
template <typename matrix_t>
inline
void
MatAdd(matrix_t& mOut, const matrix_t& m1, const matrix_t& m2)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
		{
			mOut(i, j) = m1(i, j) - m2(i, j);
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
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
	assert(m.row_size()==m.col_size() && "ERROR in Transpose(matrix_t& m): Square Matrix needed");

	typedef typename matrix_t::size_type size_type;
	matrix_t _temp;
	for(size_type i = 1; i < m.row_size(); ++i)
		for(size_type j = 0; j < i; ++j)
			_temp(i, j) = m(i, j);

	for(size_type i = 1; i < m.row_size(); ++i)
		for(size_type j = 0; j < i; ++j)
			m(i, j) = m(j, i);

	for(size_type i = 0; i < m.row_size()-1; ++i)
		for(size_type j = i+1; j < m.col_size(); ++j)
			m(i, j) = _temp(j, i);
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

/// Set each matrix entry to a scalar (componentwise)
template <typename matrix_t>
inline
void
MatSet(matrix_t& mInOut, typename matrix_t::value_type s)
{
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < mInOut.row_size(); ++i)
		for(size_type j = 0; j < mInOut.col_size(); ++j)
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
	for(size_type i = 0; i < mInOut.row_size(); ++i)
		for(size_type j = 0; j < mInOut.col_size(); ++j)
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
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
	for(size_type i = 0; i < mOut.row_size(); ++i)
		for(size_type j = 0; j < mOut.col_size(); ++j)
		{
			mOut(i, j) = m(i,j) * s;
		}
}

template <typename matrix_t>
inline
typename matrix_t::value_type
MatFrobeniusNormSq(matrix_t& m)
{
	typename matrix_t::value_type norm = 0;
	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < m.row_size(); ++i)
		for(size_type j = 0; j < m.col_size(); ++j)
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
	for(size_type j = 0; j < m.col_size(); ++j)
	{
		sum = 0;
		for(size_type i = 0; i < m.row_size(); ++i)
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
	for(size_type i = 0; i < m.row_size(); ++i)
	{
		sum = 0;
		for(size_type j = 0; j < m.col_size(); ++j)
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
	for(size_type i = 0; i < m.row_size(); ++i)
		for(size_type j = 0; j < m.col_size(); ++j)
		{
			max = (m(i,j) > max) ? m(i,j) : max;
		}
		
	return max;
}

} // end of namespace: lgmath

#endif /* __H__LGMATH__LGMATH_MATRIX_FUNCTIONS_COMMON_IMPL__ */

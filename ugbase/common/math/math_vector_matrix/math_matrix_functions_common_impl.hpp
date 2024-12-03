/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 ��7):
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
 * "Vogel, A., Reiter, S., Rupp, M., N��gel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
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

template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiplyMBT(MathMatrix<N, M, T>& mOut, const MathMatrix<N, L, T>& m1,
        const MathMatrix<M, L, T>& m2)
{
	for(size_t i = 0; i < N; ++i)
	{
		for(size_t j = 0; j < M; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < L; ++k)
			{
				mOut(i,j) += m1(i,k) * m2(j,k);
			}
		}
	}
}

template <size_t N, size_t M, size_t L, typename T>
inline void
MatMultiplyMTB(MathMatrix<N, M, T>& mOut, const MathMatrix<L, N, T>& m1,
        const MathMatrix<L, M, T>& m2)
{
	for(size_t i = 0; i < N; ++i)
	{
		for(size_t j = 0; j < M; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < L; ++k)
			{
				mOut(i,j) += m1(k,i) * m2(k,j);
			}
		}
	}
}

template <size_t N, size_t M, typename T>
inline void
MatMultiplyMBMT(MathMatrix<N, N, T>& mOut, const MathMatrix<N, M, T>& m1,
            const MathMatrix<M, M, T>& m2)
{
	MathMatrix<M, N, T> help;

	for(size_t i = 0; i < N; ++i)
		for(size_t j = 0; j < N; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < M; ++k)
			{
				help(k,j) = 0;
				for(size_t l = 0; l < M; ++l)
				{
					help(k,j) += m2(k,l) * m1(j,l);
				}

				mOut(i,j) += m1(i,k) * help(k,j);
			}
		}
}

template <size_t N, size_t M, typename T>
inline void
MatMultiplyMTBM(MathMatrix<N, N, T>& mOut, const MathMatrix<M, N, T>& m1,
            const MathMatrix<M, M, T>& m2)
{
	MathMatrix<M, N, T> help;

	for(size_t i = 0; i < N; ++i)
		for(size_t j = 0; j < N; ++j)
		{
			mOut(i,j) = 0;
			for(size_t k = 0; k < M; ++k)
			{
				help(k,j) = 0;
				for(size_t l = 0; l < M; ++l)
				{
					help(k,j) += m2(k,l) * m1(l,j);
				}

				mOut(i,j) += m1(k,i) * help(k,j);
			}
		}
}

////////////////////////////////////////////////////////////////////////////////
// "Contraction" for Matrices (note: contraction is usually known regarding tensors!)
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
			norm += m1(i,j) * m2(i,j);
		}

	return norm;
}


////////////////////////////////////////////////////////////////////////////////
// "Deviator" and trace for Matrices
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t>
inline typename matrix_t::value_type
MatDeviatorTrace(const matrix_t& m, matrix_t& dev)
{
	typename matrix_t::value_type trace = Trace(m);

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < m.num_rows(); ++i)
	{
		for(size_type j = 0; j < m.num_cols(); ++j)
		{
			dev(i,j) = m(i,j);
		}
		dev(i,i) -= 1.0 / 3.0 * trace;
	}
	return trace;
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
inline typename MathMatrix<0,0,T>::value_type
Determinant(const MathMatrix<0,0,T>& m)
{
	return 0;
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
	return	m(0,0)*m(1,1)*m(2,2)
			+ m(0,1)*m(1,2)*m(2,0)
			+ m(0,2)*m(1,0)*m(2,1)
			- m(0,0)*m(1,2)*m(2,1)
			- m(0,1)*m(1,0)*m(2,2)
			- m(0,2)*m(1,1)*m(2,0);
}

template <typename T>
inline typename MathMatrix<4,4,T>::value_type
Determinant(const MathMatrix<4,4,T>& m)
{
	return m(0,3)*m(1,2)*m(2,1)*m(3,0)-m(0,2)*m(1,3)*m(2,1)*m(3,0)
	       - m(0,3)*m(1,1)*m(2,2)*m(3,0)+m(0,1)*m(1,3)*m(2,2)*m(3,0)
           + m(0,2)*m(1,1)*m(2,3)*m(3,0)-m(0,1)*m(1,2)*m(2,3)*m(3,0)
           - m(0,3)*m(1,2)*m(2,0)*m(3,1)+m(0,2)*m(1,3)*m(2,0)*m(3,1)
           + m(0,3)*m(1,0)*m(2,2)*m(3,1)-m(0,0)*m(1,3)*m(2,2)*m(3,1)
           - m(0,2)*m(1,0)*m(2,3)*m(3,1)+m(0,0)*m(1,2)*m(2,3)*m(3,1)
           + m(0,3)*m(1,1)*m(2,0)*m(3,2)-m(0,1)*m(1,3)*m(2,0)*m(3,2)
           - m(0,3)*m(1,0)*m(2,1)*m(3,2)+m(0,0)*m(1,3)*m(2,1)*m(3,2)
           + m(0,1)*m(1,0)*m(2,3)*m(3,2)-m(0,0)*m(1,1)*m(2,3)*m(3,2)
           - m(0,2)*m(1,1)*m(2,0)*m(3,3)+m(0,1)*m(1,2)*m(2,0)*m(3,3)
           + m(0,2)*m(1,0)*m(2,1)*m(3,3)-m(0,0)*m(1,2)*m(2,1)*m(3,3)
           - m(0,1)*m(1,0)*m(2,2)*m(3,3)+m(0,0)*m(1,1)*m(2,2)*m(3,3);
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
inline typename MathMatrix<0,0,T>::value_type
GramDeterminant(const MathMatrix<0,0,T>& m)
{
	return 0;
}

template <size_t N, typename T>
inline typename MathMatrix<N,0,T>::value_type
GramDeterminant(const MathMatrix<N,0,T>& m)
{
	return 0;
}

template <size_t M, typename T>
inline typename MathMatrix<0,M,T>::value_type
GramDeterminant(const MathMatrix<0,M,T>& m)
{
	return 0;
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
inline typename MathMatrix<0,0,T>::value_type
SqrtGramDeterminant(const MathMatrix<0,0,T>& m)
{
	return 0;
}

template <size_t N, typename T>
inline typename MathMatrix<N,0,T>::value_type
SqrtGramDeterminant(const MathMatrix<N,0,T>& m)
{
	return 0;
}

template <size_t M, typename T>
inline typename MathMatrix<0,M,T>::value_type
SqrtGramDeterminant(const MathMatrix<0,M,T>& m)
{
	return 0;
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
	UG_THROW("Inverse for matrix of size "<<M<<"x"<<N<<" not implemented. You could use GeneralizedInverse for pseudo-Inverse.");
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
RightInverse(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m)
{
	//UG_STATIC_ASSERT(M <= N, pseudo_inverse_does_not_exist);
	if (M > N) // note: this is a 'static condition', it should be eliminated by the optimizer
		UG_THROW ("RightInverse: Type mismatch, cannot right-invert a MxN-matrix with M > N!");

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
RightInverse(MathMatrix<1,1>& mOut, const MathMatrix<1,1>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
RightInverse(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
RightInverse(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m){return fabs(Inverse(mOut, m));}

////////////////////////////////////////////////////////////////////////////////
// Left-Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////

template <size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
LeftInverse(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m)
{
	//UG_STATIC_ASSERT(N <= M, pseudo_inverse_does_not_exist);
	if (N > M) // note: this is a 'static condition', it should be eliminated by the optimizer
		UG_THROW ("LeftInverse: Type mismatch, cannot right-invert a MxN-matrix with M < N!");

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
LeftInverse(MathMatrix<1,1>& mOut, const MathMatrix<1,1>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<2,2,T>::value_type
LeftInverse(MathMatrix<2,2>& mOut, const MathMatrix<2,2>& m){return fabs(Inverse(mOut, m));}

template <typename T>
inline typename MathMatrix<3,3,T>::value_type
LeftInverse(MathMatrix<3,3>& mOut, const MathMatrix<3,3>& m){return fabs(Inverse(mOut, m));}

////////////////////////////////////////////////////////////////////////////////
// Generalized-Inverse of Matrix
////////////////////////////////////////////////////////////////////////////////
template<size_t N, size_t M, typename T>
inline typename MathMatrix<N,M,T>::value_type
GeneralizedInverse(MathMatrix<N,M,T>& mOut, const MathMatrix<M,N,T>& m)
{
	if(M<N){//UG_LOG("Right Inverse for matrix of size "<<M<<"x"<<N<<".");
		return RightInverse(mOut,m);
	}

	if(M>N){//UG_LOG("Left Inverse for matrix of size "<<M<<"x"<<N<<".");
		return LeftInverse(mOut,m);
	}	
	return Inverse(mOut,m);
}

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
MatDivide(matrix_t& mOut, const matrix_t& m, typename matrix_t::value_type s)
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
///	Creates a householder matrix given the orthogonal vector to the
///	householder-hypersphere through the origin.
////////////////////////////////////////////////////////////////////////////////

template <typename matrix_t, typename vector_t>
inline void
MatHouseholder(matrix_t& mOut, const vector_t& orthoVec)
{
	assert(vector_t::Size == matrix_t::RowSize);
	assert(vector_t::Size == matrix_t::ColSize);

	typename vector_t::value_type scalarProd = VecDot(orthoVec, orthoVec);

	typedef typename matrix_t::size_type size_type_mat;
	for(size_type_mat i = 0; i < mOut.num_rows(); ++i)
	{
		for(size_type_mat j = 0; j < mOut.num_cols(); ++j){
			mOut(i,j) = - 2.0/scalarProd * orthoVec[i] * orthoVec[j];
		}
		mOut(i,i) += 1.0;
	}

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

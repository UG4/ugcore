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

#ifndef __H__LGMATH__MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__
#define __H__LGMATH__MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__

#include <cmath>
#include "math_matrix.h"
#include "math_vector.h"

namespace ug
{

/// Matrix - Vector Multiplication
// vOut = m * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
MatVecMult(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v)
{
	assert(vector_t_out::Size == matrix_t::RowSize);
	assert(vector_t_in::Size == matrix_t::ColSize);

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = 0.0;
		for(size_type j = 0; j < v.size(); ++j)
		{
			vOut[i] += m(i,j) * v[j];
		}
	}
}

/// Matrix - Vector Multiplication adding to a second matrix
// vOut += m * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
MatVecMultAppend(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v)
{
	assert(vector_t_out::Size == matrix_t::RowSize);
	assert(vector_t_in::Size == matrix_t::ColSize);

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		for(size_type j = 0; j < v.size(); ++j)
		{
			vOut[i] += m(i,j) * v[j];
		}
	}
}

/// Matrix - Vector Multiplication added scaled to a second vector
// vOut += s * m * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
MatVecScaleMultAppend(vector_t_out& vOut, typename vector_t_out::value_type s, const matrix_t& m, const vector_t_in& v)
{
	assert(vector_t_out::Size == matrix_t::RowSize);
	assert(vector_t_in::Size == matrix_t::ColSize);

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		for(size_type j = 0; j < v.size(); ++j)
		{
			vOut[i] += s * m(i,j) * v[j];
		}
	}
}


/// Transposed Matrix - Vector Muliplication
// vOut = Transpose(m) * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
TransposedMatVecMult(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v)
{
	assert(vector_t_out::Size == matrix_t::ColSize);
	assert(vector_t_in::Size == matrix_t::RowSize);

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		vOut[i] = 0.0;
		for(size_type j = 0; j < v.size(); ++j)
		{
			vOut[i] += m(j,i) * v[j];
		}
	}
}

/// Transposed Matrix - Vector Muliplication
// vOut += Transpose(m) * v
template <typename vector_t_out, typename matrix_t, typename vector_t_in>
inline
void
TransposedMatVecMultAdd(vector_t_out& vOut, const matrix_t& m, const vector_t_in& v)
{
	assert(vector_t_out::Size == matrix_t::ColSize);
	assert(vector_t_in::Size == matrix_t::RowSize);

	typedef typename matrix_t::size_type size_type;
	for(size_type i = 0; i < vOut.size(); ++i)
	{
		for(size_type j = 0; j < v.size(); ++j)
		{
			vOut[i] += m(j,i) * v[j];
		}
	}
}

/**
 * Multiplication of a vector v by the Givens rotation transforming a given
 * matrix A to an upper triangular form R. After the call A stores the upper
 * triangular form R. Thus, this function computes the QR-decomposition
 * \f$ A = Q R \f$ with \f$ Q^T = Q^{-1} \f$ and performs the multiplication
 * \f$ v \gets Q^{-1} v \f$. Note that A is not necessarily a square matrix,
 * but it should have more lines than columns.
 */
template <typename matrix_t, typename vector_t>
inline
void
GivensMatVecMult (matrix_t& A, vector_t& v)
{
	typedef typename matrix_t::size_type size_type;
	typedef typename matrix_t::value_type value_type;
	
	assert (vector_t::Size == matrix_t::RowSize);
	assert (matrix_t::RowSize >= matrix_t::ColSize);
	
	value_type s, c;
	value_type d, x, y;

	for (size_type i = 0; i < matrix_t::RowSize - 1; i++) // the 1st index of the elementar rotation
		for (size_type k = matrix_t::RowSize - 1; k > i; k--) // the 2nd index of this rotation
		{
			d = A (i, i); x = A (k, i);
			
		// Computation of the entries of the elementar transformation:
			if (std::abs (d) > std::abs (x))
			{
				s = x / d;
				c = 1.0 / std::sqrt (1 + s * s);
				s *= c;
			}
			else if (x != 0) // to ensure that A (k, i) != 0; note that abs(x) >= abs(d) here!
			{
				c = d / x;
				s = 1.0 / std::sqrt (1 + c * c);
				c *= s;
			}
			else continue; // nothing to eliminate

		// Multiplication of A by the elementar transformation:
			for (size_type j = i; j < matrix_t::ColSize; j++)
			{
				x = A (i, j); y = A (k, j);
				A (i, j) = c * x + s * y;
				A (k, j) = c * y - s * x;
			}

		// Multiplication of v by the elementar transformation:
		   x = v [i]; y = v [k];
		   v [i] = c * x + s * y;
		   v [k] = c * y - s * x;
		}
}

/**
 * Multiplication of a vector v by \f$ A^{-1} \f$ fr a given matrix A
 * (\f$ v := A^{-1} v \f$) using the QR-decomposition based on the Givens
 * rotations. Note that A is not necessarily a square matrix, but it should
 * have more lines than columns. In the latter case, the result is stored
 * in v[0]...v[A.num_cols() - 1], whereas the Euclidean norm of the rest
 * is the Euclidean distance between the original v and its projection to
 * the space spanned by the columns of the matrix. If the matrix is singular
 * (i.e. its columns are linearly dependent) then an exception is thrown.
 *
 * This function can be used in the least square method and for the orthogonal
 * projection of v to the space spanned by the columns of (the original) A.
 *
 * Remark: After the call, A stores the upper triangular form of the
 * QR-decomposition, i.e. the original matrix is destroyed.
 */
template <typename matrix_t, typename vector_t>
inline
void
InvMatVecMult_byGivens (matrix_t& A, vector_t& v)
{
	typedef typename matrix_t::size_type size_type;
	typedef typename matrix_t::value_type value_type;
	
// I. Multiply 'this' by the Givens rotation:
	GivensMatVecMult (A, v);
  
// II. Computation of the result:
	size_type i = matrix_t::ColSize; // <= matrix_t::RowSize, i.e. we invert only the square block
	do
	{
		i--;
		for (size_type j = i + 1; j < matrix_t::ColSize; j++)
			v [i] -= A (i, j) * v [j];
		if (std::abs (A (i, i)) < SMALL * std::abs (v [i]))
			UG_THROW ("InvMatVecMult_byGivens: Inverting singular matrix by the Givens rotations");
		v [i] /= A (i, i);
	}
	while (i != 0);
}

/**
 * Orthogonal projection of a given vector v to the space spanned by the columns
 * of a given matrix A. The projection is written to v.
 */
template <typename matrix_t, typename vector_t>
inline
void
OrthogProjectVec (vector_t& v, const matrix_t& A)
{
	typedef typename matrix_t::size_type size_type;
	typedef typename matrix_t::value_type value_type;
	
//	I. Solve the least square problem:
	matrix_t M = A; // we do not work with the original matrix; otherwise it would be destroyed
	InvMatVecMult_byGivens (M, v);
	
//	II. Compute the linear combination of the columns:
	MathVector<matrix_t::ColSize, value_type> coeff;
	for (size_type i = 0; i < matrix_t::ColSize; i++) coeff [i] = v [i];
	MatVecMult (v, A, coeff);
}

} // end of namespace ug

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__ */


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



}// end of namespace: lgmath

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__ */

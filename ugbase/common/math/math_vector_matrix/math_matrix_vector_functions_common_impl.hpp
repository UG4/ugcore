/*
 * lgmath_matrix_vector_functions_common_impl.hpp
 *
 *  Created on: 07.07.2009
 *      Author: andreasvogel
 */

#ifndef __H__LGMATH__MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__
#define __H__LGMATH__MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__

#include <cmath>
#include "math_matrix.h"
#include "math_vector.h"

namespace ug
{

/// Matrix - Vector Muliplication
// vOut = m * v
template <int N, int M>
inline
void
MatVecMult(MathVector<N>& vOut, const MathMatrix<N,M>& m, const MathVector<M>& v)
{
	for(uint i = 0; i < MathVector<N>::Size; ++i)
	{
		vOut.coord(i) = 0.0;
		for(uint j = 0; j < MathVector<M>::Size; ++j)
		{
			vOut.coord(i) += m.entry(i,j) * v.coord(j);
		}
	}
}

/// Transposed Matrix - Vector Muliplication
// vOut = Transpose(m) * v
template <int N, int M>
inline
void
TransposedMatVecMult(MathVector<N>& vOut, const MathMatrix<M,N>& m, const MathVector<M>& v)
{
	for(uint i = 0; i < MathVector<N>::Size; ++i)
	{
		vOut.coord(i) = 0.0;
		for(uint j = 0; j < MathVector<M>::Size; ++j)
		{
			vOut.coord(i) += m.entry(j,i) * v.coord(j);
		}
	}
}



}// end of namespace: lgmath

#endif /* __H__LGMATH__LGMATH_MATRIX_VECTOR_FUNCTIONS_COMMON_IMPL__ */

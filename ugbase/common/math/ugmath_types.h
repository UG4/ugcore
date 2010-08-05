//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d12

#ifndef __H__UGMATH_TYPES__
#define __H__UGMATH_TYPES__

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	This header defines common matrix and vector types

////////////////////////////////////////////////////////////////////////
//	include vector and matrix classes
#include "math_vector_matrix/math_vector.h"
#include "math_vector_matrix/math_matrix.h"
#include "math_vector_matrix/math_tensor.h"

////////////////////////////////////////////////////////////////////////
//	useful typedefs
namespace ug
{
typedef MathVector<2> vector2;
typedef MathVector<3> vector3;
typedef MathVector<4> vector4;

typedef MathMatrix<2,2> matrix22;
typedef MathMatrix<3,3> matrix33;
typedef MathMatrix<4,4> matrix44;
}

#endif

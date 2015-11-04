#ifndef __H__UGMATH_TYPES__
#define __H__UGMATH_TYPES__

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	This header defines common matrix and vector types

////////////////////////////////////////////////////////////////////////
//	include vector and matrix classes
#include "math_vector_matrix/math_vector.h"
#include "math_vector_matrix/math_matrix.h"
#include "math_vector_matrix/math_symmetric_matrix.h"
#include "math_vector_matrix/math_tensor.h"

////////////////////////////////////////////////////////////////////////
//	useful typedefs
namespace ug
{

/// \addtogroup vectors
/// \{

/// a 1d vector
typedef MathVector<1, number> vector1;

/// a 2d vector
typedef MathVector<2, number> vector2;

/// a 3d vector
typedef MathVector<3, number> vector3;

/// a 4d vector
typedef MathVector<4, number> vector4;

/// a 2x2 matrix
typedef MathMatrix<2,2, number> matrix22;

/// a 3x3 matrix
typedef MathMatrix<3,3, number> matrix33;

/// a 4x4 matrix
typedef MathMatrix<4,4, number> matrix44;

// end group vectors
/// \}

}

#endif

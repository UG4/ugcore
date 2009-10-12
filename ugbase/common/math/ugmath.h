//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d12

#ifndef __H__UGMATH__
#define __H__UGMATH__

#include <cmath>
////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	This header defines common math and defines some vector-types.

////////////////////////////////////////////////////////////////////////
//	includes of useful mathematical methods
#include "misc/math_util.h"


////////////////////////////////////////////////////////////////////////
//	include vector and matrix classes
#include "math_vector_matrix/math_vector.h"
#include "math_vector_matrix/math_matrix.h"

////////////////////////////////////////////////////////////////////////
// Functions
#include "math_vector_matrix/math_vector_functions.h"
#include "math_vector_matrix/math_matrix_functions.h"
#include "math_vector_matrix/math_matrix_vector_functions.h"

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

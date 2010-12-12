/**
 * \file ugmath_types.h
 * \brief typedefs for ugmath
 *
 * \author Sebastian Reiter
 * \date  y08 m11 d12
 *
 * This file contains typedefs for ug math.
 */

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

/**
 * Abbreviations of small vectors
 *
 * \defgroup vectors Vectors
 */

/// \addtogroup vectors
///@{

/// a 1d vector
typedef MathVector<1> vector1;

/// a 2d vector
typedef MathVector<2> vector2;

/**
 * A 3d vector
 */
typedef MathVector<3> vector3;

/// a 4d vector
typedef MathVector<4> vector4;

/// a 2x2 matrix
typedef MathMatrix<2,2> matrix22;

/**
 * A 3x3 matrix
 */
typedef MathMatrix<3,3> matrix33;

/// a 4x4 matrix
typedef MathMatrix<4,4> matrix44;

///@}

}

#endif

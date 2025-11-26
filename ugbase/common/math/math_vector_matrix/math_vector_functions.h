/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__COMMON__MATH_VECTOR_FUNCTIONS__
#define __H__COMMON__MATH_VECTOR_FUNCTIONS__

#include "math_vector.h"

namespace ug
{

/// \addtogroup vectors
/// \{

////////////////////////////////////////////////////////////////
///	Copy contents between vectors of possibly different types.
/**	If the target vector is bigger than the source vector, additional entries
 * will be set to the scalar given in the fill argument. If it is smaller,
 * then the additional entries in the source vector are ignored.*/
template <typename vector_target_t, typename vector_source_t>
void VecCopy(vector_target_t& target, const vector_source_t& source,
			 typename vector_target_t::value_type fill);

////////////////////////////////////////////////////////////////
// Append Vector

///	adds a MathVector<N> to a second one
// vOut += v1
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1);

///	adds two MathVector<N>s and adds the result to a third one
// vOut = v1 + v2
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1, const vector_t& v2);

///	adds two MathVector<N>s and adds the result to a fourth one
// vOut += v1 + v2 + v3
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1, const vector_t& v2,
							const vector_t& v3);

///	adds two MathVector<N>s and adds the result to a fifth one
// vOut += v1 + v2 + v3 + v4
template <typename vector_t>
inline
void
VecAppend(vector_t& vOut, const vector_t& v1, const vector_t& v2,
							const vector_t& v3, const vector_t& v4);

////////////////////////////////////////////////////////////////
// Scale and Append Vectors

/// Scales a Vector and adds it to a second vector
// vOut += s1*v1
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1);

/// Scales two Vectors, adds them and adds the sum to a third vector
// vOut += s1*v1 + s2*v2
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2);

/// Scales three Vectors, adds them and adds the sum to a fourth vector
// vOut += s1*v1 + s2*v2 + s3*v3
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3);

/// Scales four Vectors, adds them and adds the sum to a fifth vector
// vOut += s1*v1 + s2*v2 + s3*v3 + s4*v4
template <typename vector_t>
inline
void
VecScaleAppend(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3,
								 typename vector_t::value_type s4, const vector_t& v4);


////////////////////////////////////////////////////////////////
// Addition of Vectors

///	adds two MathVector<N>s and stores the result in a third one
// vOut = v1 + v2
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2);

///	adds three MathVector<N>s and stores the result in a fourth one
// vOut = v1 + v2 + v3
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2,
							const vector_t& v3);

///	adds four MathVector<N>s and stores the result in a firth one
// vOut = v1 + v2 + v3 + v4
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v1, const vector_t& v2,
							const vector_t& v3, const vector_t& v4);

////////////////////////////////////////////////////////////////
// Subtraction of Vectors

///	subtracts v2 from v1 and stores the result in a vOut
// vOut = v1 - v2
template <typename vector_t>
inline
void
VecSubtract(vector_t& vOut, const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////
// Component-wise exponentiation of a vector
template <typename vector_t>
inline
void
VecPow(vector_t& vOut, const vector_t& v1, typename vector_t::value_type s);


////////////////////////////////////////////////////////////////
// Scaling of a Vector

///	scales a MathVector<N>
// vOut = s * v
template <typename vector_t>
inline
void
VecScale(vector_t& vOut, const vector_t& v, typename vector_t::value_type s);

////////////////////////////////////////////////////////////////
// Scaled Addition of Vectors

/// Scales two Vectors, adds them and returns the sum in a third vector
// vOut = s1*v1 + s2*v2
template <typename vector_t>
inline
void
VecScaleAdd(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2);

/// Scales three Vectors, adds them and returns the sum in a fourth vector
// vOut = s1*v1 + s2*v2 + s3*v3
template <typename vector_t>
inline
void
VecScaleAdd(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3);

/// Scales four Vectors, adds them and returns the sum in a fifth vector
// vOut = s1*v1 + s2*v2 + s3*v3 + s4*v4
template <typename vector_t>
inline
void
VecScaleAdd(vector_t& vOut, typename vector_t::value_type s1, const vector_t& v1,
								 typename vector_t::value_type s2, const vector_t& v2,
								 typename vector_t::value_type s3, const vector_t& v3,
								 typename vector_t::value_type s4, const vector_t& v4);

////////////////////////////////////////////////////////////////
// Interpolation of Vectors

//	performs linear interpolation between two Vectors.
template <typename vector_t>
inline
void
VecInterpolateLinear(vector_t& vOut, const vector_t& v1, const vector_t& v2,
						typename vector_t::value_type interpAmount);

////////////////////////////////////////////////////////////////
// Length of Vectors

///	returns the squared length of v. Faster than VecLength.
template <typename vector_t>
inline
typename vector_t::value_type
VecLengthSq(const vector_t& v);

///	returns the length of v. Slower than VecLengthSq.
template <typename vector_t>
inline
typename vector_t::value_type
VecLength(const vector_t& v);

////////////////////////////////////////////////////////////////
// Distance of Vectors

///	returns the squared distance of two vector_ts.
template <typename vector_t>
inline
typename vector_t::value_type
VecDistanceSq(const vector_t& v1, const vector_t& v2);

///	returns the distance of two vector_ts.
template <typename vector_t>
inline
typename vector_t::value_type
VecDistance(const vector_t& v1, const vector_t& v2);


////////////////////////////////////////////////////////////////
// Dot Product
///	returns the dot-product of two vector_ts
template <typename vector_t>
inline
typename vector_t::value_type
VecDot(const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////
// Angle
///	returns the angle between two vectors in radiants
template <typename vector_t>
inline
typename vector_t::value_type
VecAngle(const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////
// Angle
///	returns the angle between two vectors of length 1 in radiants
/** This method assumes that v1 and v2 are normalized vectors (both have length 1).
 * No checks are performed!*/
template <typename vector_t>
inline
typename vector_t::value_type
VecAngleNorm(const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////
// Cross Product
///	calculates the cross product of two Vectors of dimension 3. It makes no sense to use VecCross for vector_ts that have not dimension 3.
template <typename vector_t>
inline
void
VecCross(vector_t& vOut, const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////
// Generalized cross Product
///	calculates the usual cross product in 3d, and the (det, 0) vector as a cross product in 2d
template <size_t dim>
inline void GenVecCross
(
	MathVector<dim> & result, ///< = v_1 X v_2
	const MathVector<dim> & v_1, const MathVector<dim> & v_2
);

////////////////////////////////////////////////////////////////
// Normalize a Vector
///	scales a vector_t to unit length
template <typename vector_t>
inline
void
VecNormalize(vector_t& vOut, const vector_t& v);

////////////////////////////////////////////////////////////////
// CalculateTriangleNormalNoNormalize
///	Calculates a triangle-normal in 3d (no normalization is performed).
//TODO: This method should not be part of the core-math methods.
//		Instead it should be moved to a seperate math-util file.
template <typename vector_t>
inline
void
CalculateTriangleNormalNoNormalize(vector_t& vOut, const vector_t& v1,
						const vector_t& v2, const vector_t& v3);

////////////////////////////////////////////////////////////////
// CalculateTriangleNormal
///	Calculates a triangle-normal in 3d (output has length 1).
//TODO: This method should not be part of the core-math methods.
//		Instead it should be moved to a seperate math-util file.
template <typename vector_t>
inline
void
CalculateTriangleNormal(vector_t& vOut, const vector_t& v1,
						const vector_t& v2, const vector_t& v3);


////////////////////////////////////////////////////////////////
// Operations with scalar

/// Set each vector component to scalar (componentwise)
template <typename vector_t>
inline
void
VecSet(vector_t& vInOut, typename vector_t::value_type s);

/// Add a scalar to a vector (componentwise)
template <typename vector_t>
inline
void
VecAdd(vector_t& vOut, const vector_t& v, typename vector_t::value_type s);

/// Subtract a scalar from a vector (componentwise)
template <typename vector_t>
inline
void
VecSubtract(vector_t& vOut, const vector_t& v, typename vector_t::value_type s);

////////////////////////////////////////////////////////////////
// Norms
template <typename vector_t>
inline
typename vector_t::value_type
VecTwoNorm(const vector_t& v);

template <typename vector_t>
inline
typename vector_t::value_type
VecTwoNormSq(const vector_t& v);

template <typename vector_t>
inline
typename vector_t::value_type
VecOneNorm(const vector_t& v);

template <typename vector_t>
inline
typename vector_t::value_type
VecPNorm(const vector_t& v, unsigned int p);

template <typename vector_t>
inline
typename vector_t::value_type
VecMaxNorm(const vector_t& v);

template <typename vector_t>
inline
typename vector_t::value_type
VecInftyNorm(const vector_t& v);

// end group vectors
/// \}


////////////////////////////////////////////////////////////////
// Elementwise operations

/// component-wise product: vOut_i = v1_i*v2_i
template <typename vector_t>
inline
void
VecElemProd(vector_t& vOut, const vector_t& v1, const vector_t& v2);

///	component-wise square root: \f$o_i=\sqrt {v_i}\f$
template <typename vector_t>
inline
void
VecElemSqrt(vector_t& vOut, const vector_t& v1);

///	component-wise comparison of two vectors (in the absolute values)
template <typename vector_t>
inline
bool
VecAbsIsLess(const vector_t& v1, const vector_t& v2);

///	component-wise comparison of a vector (in the absolute values) with a given number
template <typename vector_t>
inline
bool
VecAbsIsLess(const vector_t& v1, const typename vector_t::value_type s);

/// checks if the given point is in the bounding box given by two other points
template <typename vector_t>
inline
bool
VecIsInBB(const vector_t& v, const vector_t& low, const vector_t& high);

}//	end of namespace

////////////////////////////////////////////////////////////////////////
//	include a general, but not very fast implementation of the declared methods above.
#include "math_vector_functions_common_impl.hpp"

#endif

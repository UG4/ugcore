//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d14

#ifndef __H__UGMATH__MATH_UTIL__
#define __H__UGMATH__MATH_UTIL__

#include <cstdlib>
#include "../math_vector_matrix/math_vector.h"
#include "../math_vector_matrix/math_vector_functions.h"

namespace ug
{

const number SMALL = 1.0e-12;

////////////////////////////////////////////////////////////////////////
//	urand
///	uniform distributed random numbers in [lowerBound, upperBound[. Use srand to set a seed.
template <class TNumber>
TNumber
urand(TNumber lowerBound, TNumber upperBound);

////////////////////////////////////////////////////////////////////////
//	clip
///	clips a number to the given interval [lowerBound, upperBound].
template <class TNumber>
TNumber
clip(TNumber val, TNumber lowerBound, TNumber upperBound);



////////////////////////////////////////////////////////////////////////
///	finds the projection of v onto the line defined by v0 and v1
/**
 * projects v onto the line defined by v0 and v1.
 * The projected point is returned in vOut.
 *
 * returns s so that vOut = (1.0-s)*v0 + s*v1
 */
template <class vector_t>
number DropAPerpendicular(vector_t& vOut, const vector_t& v0,
							const vector_t& v1, const vector_t& v);

////////////////////////////////////////////////////////////////////////
///	projects v onto the plane defined by the point p and the planes normal n.
/**	The result is written to vOut.*/
template <class vector_t>
void ProjectPointToPlane(vector_t& vOut, const vector_t& v,
						const vector_t& p, const vector_t& n);

////////////////////////////////////////////////////////////////////////
///	calculates the are of the triangle defined by p1, p2 and p3
template <class vector_t>
number TriangleArea(const vector_t& p1, const vector_t& p2,
							 const vector_t& p3);

////////////////////////////////////////////////////////////////////////
///	the returned degree lies between 0 and 1. The closer to 1 the better.
/**
 * returns the minimal dot-product of each normal of the triangle-corners
 * with the triangles normal.
 */
template <class vector_t>
number GeometricApproximationDegree(vector_t& n1, vector_t& n2, vector_t& n3,
									vector_t& tn);

////////////////////////////////////////////////////////////////////////
///	returns a value between 0 and 1. The higher the better.
/**	The quality is measured by comparing the triangles area with the
 *	sum of the squared length of each side.
 *
 *	a triangle whose sides have all zero length is considered to be
 *	a bad triangle - 0 is returned.*/
template <class vector_t>
number TriangleQuality_Area(const vector_t& p1, const vector_t& p2,
							const vector_t& p3);

}//	end of namespace

////////////////////////////////////////////////
//	include implementation
#include "math_util_impl.hpp"

#endif

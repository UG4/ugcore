/*
 * geometry_util.h
 *
 *  Created on: 07.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__COMMON__GEOMETRY_UTIL__
#define __H__LIB_DISCRETIZATION__COMMON__GEOMETRY_UTIL__

#include <vector>
#include <cmath>

#include "common/common.h"
#include "lib_discretization/reference_element/reference_element.h"

namespace ug{

///////////////////////////////////////////////////////////////
/// Volume of an Element in a given Dimension
template <typename TRefElem, int TWorldDim>
inline number ElementSize(const MathVector<TWorldDim>* vCornerCoords);

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// 1D Reference Elements
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
/// Volume of a Line in 1d
/**
 * This function returns the volume of a line in 1d.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (2 corners)
 * \return 		number			Volume of Line
 */
template <>
inline number ElementSize<ReferenceEdge, 1>(const MathVector<1>* vCornerCoords)
{
	return VecDistance(vCornerCoords[0], vCornerCoords[1]);
}

///////////////////////////////////////////////////////////////
/// Volume of a Line in 2d
/**
 * This function returns the volume of a line in 2d.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (2 corners)
 * \return 		number			Volume of Line
 */
template <>
inline number ElementSize<ReferenceEdge, 2>(const MathVector<2>* vCornerCoords)
{
	return VecDistance(vCornerCoords[0], vCornerCoords[1]);
}

///////////////////////////////////////////////////////////////
/// Volume of a Line in 3d
/**
 * This function returns the volume of a line in 3d.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (2 corners)
 * \return 		number			Volume of Line
 */
template <>
inline number ElementSize<ReferenceEdge, 3>(const MathVector<3>* vCornerCoords)
{
	return VecDistance(vCornerCoords[0], vCornerCoords[1]);
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// 2D Reference Elements
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
/// Volume of a Triangle in 2d
/**
 * This function returns the volume of a triangle in 2d.
 * The Volume is computed via:
 * 	F = 1/2 * fabs( (x2-x1)*(y1-y0) - (x1-x0)*(y2-y0) )
 *
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (3 corners)
 * \return 		number			Volume of Triangle
 */
template <>
inline number ElementSize<ReferenceTriangle, 2>(const MathVector<2>* vCornerCoords)
{
	return(0.5*fabs((vCornerCoords[1][1]-vCornerCoords[0][1])*(vCornerCoords[2][0]-vCornerCoords[0][0])
				   -(vCornerCoords[1][0]-vCornerCoords[0][0])*(vCornerCoords[2][1]-vCornerCoords[0][1])));
}

///////////////////////////////////////////////////////////////
/// Volume of a Triangle in 3d
/**
 * This function returns the volume of a triangle in 3d.
 * The Volume is computed via:
 * 	F = 1/2 * | (c - a) x (b - a) | (Here, a,b,c are the vectors of the corners)
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (3 corners)
 * \return 		number			Volume of Triangle
 */
template <>
inline number ElementSize<ReferenceTriangle, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x10, x20, n;

	VecSubtract(x10, vCornerCoords[1], vCornerCoords[0]);
	VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);
	VecCross(n, x10, x20);

	return(0.5 * VecTwoNorm(n) );
}

///////////////////////////////////////////////////////////////
/// Volume of a Quadrilateral in 2d
/**
 * This function returns the volume of a quadrilateral in 2d.
 * The Volume is computed via:
 * 	F = 1/2 * | (y3-y1)*(x2-x0) - (x3-x1)*(y2-y0) |
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (4 corners)
 * \return 		number			Volume of Quadrilateral
 */
template <>
inline number ElementSize<ReferenceQuadrilateral, 2>(const MathVector<2>* vCornerCoords)
{
	return( 0.5*fabs( (vCornerCoords[3][1]-vCornerCoords[1][1])*(vCornerCoords[2][0]-vCornerCoords[0][0])
					 -(vCornerCoords[3][0]-vCornerCoords[1][0])*(vCornerCoords[2][1]-vCornerCoords[0][1]) ) );
}

///////////////////////////////////////////////////////////////
/// Volume of a Quadrilateral in 3d
/**
 * This function returns the volume of a quadrilateral in 3d.
 * The Volume is computed via:
 * 	F = 1/2 * | (c - a) x (d - b) |  (Here, a,b,c,d are the vectors of the corners)
 *
 * The corner coordinates must be given in counter - clockwise order
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (4 corners)
 * \return 		number			Volume of Quadrilateral
 */
template <>
inline number ElementSize<ReferenceQuadrilateral, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x20, x31, n;

	VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);
	VecSubtract(x31, vCornerCoords[3], vCornerCoords[1]);
	VecCross(n, x20, x31);

	return(0.5 * VecTwoNorm(n) );
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// 3D Reference Elements
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
/// Volume of a Tetrahedron in 3d
/**
 * This function returns the volume of a tetrahedron in 3d.
 * The Volume is computed via:
 * 	V = 1/6 * | ( (b - a) x (c - a) ) * (d - a) |  (Here, a,b,c,d are the vectors of the corners)
 *
 * This is motivated by the formula: V = 1/3 * (S * h) with
 * - S is the area of the base
 * - h is the height
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (4 corners)
 * \return 		number			Volume of Tetrahedron
 */
template <>
inline number ElementSize<ReferenceTetrahedron, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x10, x20, x30, n;

	VecSubtract(x10, vCornerCoords[1], vCornerCoords[0]);
	VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);
	VecSubtract(x30, vCornerCoords[3], vCornerCoords[0]);
	VecCross(n, x10, x20);

	return (1./6.) * fabs ( VecDot(n, x30) );
}

///////////////////////////////////////////////////////////////
/// Volume of a Pyramid in 3d
/**
 * This function returns the volume of a pyramid in 3d.
 * The volume is computed via: V = 1/3 * (S * h) with
 * - S is the area of the base
 * - h is the height
 *
 * The corner coordinates must be given as prescribed by the reference element
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (5 corners)
 * \return 		number			Volume of Pyramid
 */
template <>
inline number ElementSize<ReferencePyramid, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x20, x31, x40, n;

	VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);
	VecSubtract(x31, vCornerCoords[3], vCornerCoords[1]);
	VecSubtract(x40, vCornerCoords[4], vCornerCoords[0]);
	VecCross(n, x20, x31);

	return (1./6.) * VecDot(n, x40);
}

///////////////////////////////////////////////////////////////
/// Volume of a Prism in 3d
/**
 * This function returns the volume of a prism in 3d. Therefore, the
 * element is devided into pyramid {x0, x_1, x_4, x_3; x_5} and  a
 * tetrahedron {x_0, x_1, x_2; x_5}, whose volumes are computed and added.
 *
 * The corner coordinates must be given as prescribed by the reference element
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (6 corners)
 * \return 		number			Volume of Prism
 */
template <>
inline number ElementSize<ReferencePrism, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x40, x13, x10, x20, x50, m, n;

	VecSubtract(x40, vCornerCoords[4], vCornerCoords[0]);
	VecSubtract(x13, vCornerCoords[1], vCornerCoords[3]);
	VecSubtract(x10, vCornerCoords[1], vCornerCoords[0]);
	VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);

	// height for both subelements
	VecSubtract(x50, vCornerCoords[5], vCornerCoords[0]);

	VecCross(m, x40, x13); // base of pyramid
	VecCross(n, x10, x20); // base of tetrahedron

	// n = n + m
	VecAppend(n, m);

	return (1./6.) * VecDot(n, x50);
}

///////////////////////////////////////////////////////////////
/// Volume of a Hexahedron in 3d
/**
 * This function returns the volume of a hexahedron in 3d. Therefore, the
 * element is devided into two prisms:
 *  1) {x0, x_1, x_2, x_4; x_5, x_6} and
 *  2) {x0, x_2, x_3, x_4; x_6, x_7}
 *
 * The corner coordinates must be given as prescribed by the reference element
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (8 corners)
 * \return 		number			Volume of Hexahedron
 */
template <>
inline number ElementSize<ReferenceHexahedron, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x50, x14, x10, x20, x60, m1, n1;
	MathVector<3>      x24,      x30, x70, m2, n2;

	// 1. prism
	VecSubtract(x50, vCornerCoords[5], vCornerCoords[0]);
	VecSubtract(x14, vCornerCoords[1], vCornerCoords[4]);
	VecSubtract(x10, vCornerCoords[1], vCornerCoords[0]);
	VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);
	VecSubtract(x60, vCornerCoords[6], vCornerCoords[0]);

	VecCross(m1, x50, x14); // base of pyramid
	VecCross(n1, x10, x20); // base of tetrahedron

	// n1 = n1 + m1
	VecAppend(n1, m1);

	// 2. prism
	//VecSubtract(x60, vCornerCoords[6], vCornerCoords[0]);
	VecSubtract(x24, vCornerCoords[2], vCornerCoords[4]);
	//VecSubtract(x20, vCornerCoords[2], vCornerCoords[0]);
	VecSubtract(x30, vCornerCoords[3], vCornerCoords[0]);
	VecSubtract(x70, vCornerCoords[7], vCornerCoords[0]);

	VecCross(m2, x60, x24); // base of pyramid
	VecCross(n2, x20, x30); // base of tetrahedron

	// n2 = n2 + m2
	VecAppend(n2, m2);

	return (1./6.) * (VecDot(n1, x60) + VecDot(n2, x70));
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// Normals on Elements
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////
/// Normal to an Element in a given Dimension
template <typename TRefElem, int TWorldDim>
inline void ElementNormal(MathVector<TWorldDim>& normalOut, const MathVector<TWorldDim>* vCornerCoords);

///////////////////////////////////////////////////////////////
/// Normal to a Point in 1d
/**
 * This function returns the normal of a point in 1d.
 * This always 1.0.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceVertex, 1>(MathVector<1>& normalOut, const MathVector<1>* vCornerCoords)
{
	normalOut[0] = 1.0;
}

///////////////////////////////////////////////////////////////
/// Normal to a Line in 2d
/**
 * This function returns the normal of a line in 2d.
 * The normal is in clockwise direction to the vector (x1-x0).
 * The norm of the normal is the size of the line.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceEdge, 2>(MathVector<2>& normalOut, const MathVector<2>* vCornerCoords)
{
	MathVector<2> diff(vCornerCoords[1]);
	diff -= vCornerCoords[0];

	normalOut[0] = diff[1];
	normalOut[1] = -diff[0];
}

///////////////////////////////////////////////////////////////
/// Normal to a Triangle in 3d
/**
 * This function returns the normal of a triangle in 3d.
 * The orientation is right-handed (i.e. forming with the four fingers of the
 * right hand a circle in direction of the corner nodes, the thumb points in
 * the direction of the normal)
 *
 * The norm of the normal is the size of the triangle.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceTriangle, 3>(MathVector<3>& normalOut, const MathVector<3>* vCornerCoords)
{
	MathVector<3> a, b;
	VecSubtract(a, vCornerCoords[1], vCornerCoords[0]);
	VecSubtract(b, vCornerCoords[2], vCornerCoords[0]);
	VecCross(normalOut, a,b);
	VecScale(normalOut, normalOut, 0.5);
}

///////////////////////////////////////////////////////////////
/// Normal to a Quadrilateral in 3d
/**
 * This function returns the normal of a quadrilateral in 3d.
 * The orientation is right-handed (i.e. forming with the four fingers of the
 * right hand a circle in direction of the corner nodes, the thumb points in
 * the direction of the normal)
 *
 * The norm of the normal is the size of the quadrilateral.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceQuadrilateral, 3>(MathVector<3>& normalOut, const MathVector<3>* vCornerCoords)
{
	MathVector<3> a, b;
	VecSubtract(a, vCornerCoords[2], vCornerCoords[0]);
	VecSubtract(b, vCornerCoords[3], vCornerCoords[1]);
	VecCross(normalOut, a,b);
	VecScale(normalOut, normalOut, 0.5);
}



} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__COMMON__GEOMETRY_UTIL__ */

/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__COMMON__GEOMETRY_UTIL__
#define __H__UG__LIB_DISC__COMMON__GEOMETRY_UTIL__

#include <vector>
#include <cmath>

#include "common/common.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/reference_element/element_list_traits.h"
#include "lib_disc/domain_traits.h"
#include "common/util/provider.h"

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
 * element is divided into pyramid {x0, x_1, x_4, x_3; x_5} and  a
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
 * element is divided into two prisms:
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
/// Volume of an Octahedron in 3d
/**
 * This function returns the volume of an octhedron in 3d
 * by calculating the volumes of the upper and lower pyramid
 * the octahedron consists of.
 * The pyramidal volumes are computed via: V = 1/3 * (S * h) with
 * - S is the area of the base
 * - h is the height
 *
 * The corner coordinates must be given as prescribed by the reference element
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates (6 corners)
 * \return 		number			Volume of Octahedron
 */
template <>
inline number ElementSize<ReferenceOctahedron, 3>(const MathVector<3>* vCornerCoords)
{
	MathVector<3> x31, x42, x51, n;
	MathVector<3> 			x01;

	VecSubtract(x31, vCornerCoords[3], vCornerCoords[1]);
	VecSubtract(x42, vCornerCoords[4], vCornerCoords[2]);
	VecSubtract(x51, vCornerCoords[5], vCornerCoords[1]);
	VecCross(n, x31, x42);

	//VecSubtract(x31, vCornerCoords[3], vCornerCoords[1]);
	//VecSubtract(x42, vCornerCoords[4], vCornerCoords[2]);
	VecSubtract(x01, vCornerCoords[0], vCornerCoords[1]);
	//VecCross(n, x31, x42);

	number volTopPyr 	= (1./6.) * fabs(VecDot(n, x51));
	number volBottomPyr = (1./6.) * fabs(VecDot(n, x01));

	return volTopPyr + volBottomPyr;
}

///////////////////////////////////////////////////////////////
//	run-time size of element
///////////////////////////////////////////////////////////////

template <int dim>
inline number ElementSize(ReferenceObjectID roid, const MathVector<dim>* vCornerCoords);

template <>
inline number ElementSize<1>(ReferenceObjectID roid, const MathVector<1>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: return 1.0;
		case ROID_EDGE: return ElementSize<ReferenceEdge, 1>(vCornerCoords);
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 1.");
	}
}

template <>
inline number ElementSize<2>(ReferenceObjectID roid, const MathVector<2>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: return 1.0;
		case ROID_EDGE: return ElementSize<ReferenceEdge, 2>(vCornerCoords);
		case ROID_TRIANGLE: return ElementSize<ReferenceTriangle, 2>(vCornerCoords);
		case ROID_QUADRILATERAL: return ElementSize<ReferenceQuadrilateral, 2>(vCornerCoords);
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 2.");
	}
}

template <>
inline number ElementSize<3>(ReferenceObjectID roid, const MathVector<3>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: return 1.0;
		case ROID_EDGE: return ElementSize<ReferenceEdge, 3>(vCornerCoords);
		case ROID_TRIANGLE: return ElementSize<ReferenceTriangle, 3>(vCornerCoords);
		case ROID_QUADRILATERAL: return ElementSize<ReferenceQuadrilateral, 3>(vCornerCoords);
		case ROID_TETRAHEDRON: return ElementSize<ReferenceTetrahedron, 3>(vCornerCoords);
		case ROID_PYRAMID: return ElementSize<ReferencePyramid, 3>(vCornerCoords);
		case ROID_PRISM: return ElementSize<ReferencePrism, 3>(vCornerCoords);
		case ROID_HEXAHEDRON: return ElementSize<ReferenceHexahedron, 3>(vCornerCoords);
		case ROID_OCTAHEDRON: return ElementSize<ReferenceOctahedron, 3>(vCornerCoords);
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 3.");
	}
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
/// Normal to a Point in 2d
/**
 * This function returns the normal of a point in 2d.
 * This can only be understood as the normal on a corner of an edge in 2d.
 * Therefore, it will be assumed that vCornerCoords has exactly 2 entries
 * defining the edge. Then the normal is the one-dim. normal on the first vertex
 * embedded in the 1D subspace defined by the edge.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceVertex, 2>(MathVector<2>& normalOut, const MathVector<2>* vCornerCoords)
{
	try
	{
		normalOut = vCornerCoords[0];
		normalOut -= vCornerCoords[1];
		VecNormalize(normalOut, normalOut);
	}
	UG_CATCH_THROW("Element normal can not be computed. Maybe there are not enough vertices to work on.")
}

///////////////////////////////////////////////////////////////
/// Normal to a Point in 3d
/**
 * This function returns the normal of a point in 3d.
 * This can only be understood as the normal on a corner of an edge in 3d.
 * Therefore, it will be assumed that vCornerCoords has exactly 2 entries
 * defining the edge. Then the normal is the one-dim. normal on the first vertex
 * embedded in the 1D subspace defined by the edge.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceVertex, 3>(MathVector<3>& normalOut, const MathVector<3>* vCornerCoords)
{
	try
	{
		normalOut = vCornerCoords[0];
		normalOut -= vCornerCoords[1];
		VecNormalize(normalOut, normalOut);
	}
	UG_CATCH_THROW("Element normal can not be computed. Maybe there are not enough vertices to work on.")
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
/// Normal to a Line in 3d
/**
 * This function returns the normal of a line in 3d.
 * This can only be understood as the normal on an edge in a 2d manifold
 * defined by at least a triangle. Therefore, it is assumed that vCornerCoords
 * contains at least three vertices.
 * The normal is computed as the outer normal on the first edge (v2-v1) of the triangle
 * and embedded in the 2d subspace defined by the triangle.
 * The norm of the normal is the size of the line.
 *
 * \param[in]	vCornerCoords	Vector of corner coordinates
 * \param[out]	normalOut		Normal
 */
template <>
inline void ElementNormal<ReferenceEdge, 3>(MathVector<3>& normalOut, const MathVector<3>* vCornerCoords)
{
	try
	{
		// normal on triangle
		MathVector<3> edge0, edge1;
		VecSubtract(edge0, vCornerCoords[1], vCornerCoords[0]);
		VecSubtract(edge1, vCornerCoords[2], vCornerCoords[1]);
		VecCross(normalOut, edge0, edge1);
		// normal an edge is edge x normal on triangle
		VecCross(normalOut, edge0, normalOut);
		// scale
		VecNormalize(normalOut, normalOut);
		VecScale(normalOut, normalOut, VecTwoNorm(edge0));
	}
	UG_CATCH_THROW("Element normal can not be computed. Maybe there are not enough vertices to work on.")
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

///////////////////////////////////////////////////////////////
//	run-time normal of element
///////////////////////////////////////////////////////////////

template <int dim>
inline void ElementNormal(ReferenceObjectID roid, MathVector<dim>& normalOut, const MathVector<dim>* vCornerCoords);

template <>
inline void ElementNormal<1>(ReferenceObjectID roid, MathVector<1>& normalOut, const MathVector<1>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: ElementNormal<ReferenceVertex, 1>(normalOut, vCornerCoords); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 1.");
	}
}

template <>
inline void ElementNormal<2>(ReferenceObjectID roid, MathVector<2>& normalOut, const MathVector<2>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: ElementNormal<ReferenceVertex, 2>(normalOut, vCornerCoords); return;
		case ROID_EDGE: ElementNormal<ReferenceEdge, 2>(normalOut, vCornerCoords); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 2.");
	}
}

template <>
inline void ElementNormal<3>(ReferenceObjectID roid, MathVector<3>& normalOut, const MathVector<3>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: ElementNormal<ReferenceVertex, 3>(normalOut, vCornerCoords); return;
		case ROID_EDGE: ElementNormal<ReferenceEdge, 3>(normalOut, vCornerCoords); return;
		case ROID_TRIANGLE: ElementNormal<ReferenceTriangle, 3>(normalOut, vCornerCoords); return;
		case ROID_QUADRILATERAL: ElementNormal<ReferenceQuadrilateral, 3>(normalOut, vCornerCoords); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 3.");
	}
}

///////////////////////////////////////////////////////////////
/// Normal to a side of an Element in a given Dimension
/**
 * This function computes the outer normal to a given side
 * of an element. The Euclidean norm of the normal is the
 * area of the side. Note that the normal is computed in
 * the same dimensionality as the reference element itself.
 *
 * \param[in] side			index of the side of the element
 * \param[in] vCornerCoords	array of the global coordinates of the corners
 * \param[out] normalOut	the computed normal
 */
template <typename TRefElem, int TWorldDim>
inline void SideNormal(MathVector<TWorldDim>& normalOut, int side, const MathVector<TWorldDim>* vCornerCoords)
{
	static constexpr int dim = TRefElem::dim;
	static constexpr int maxSideCorners = element_list_traits<typename domain_traits<dim-1>::DimElemList>::maxCorners;
	
//	Get own reference element and the side roid:
	auto & rRefElem = (TRefElem&) ReferenceElementProvider::get(TRefElem::REFERENCE_OBJECT_ID);
	ReferenceObjectID sideRoid = rRefElem.roid(dim-1,side);
	
//	Get the coordinates of the vertices:
	MathVector<TWorldDim> vSideCorner [(dim == TWorldDim)? maxSideCorners : maxSideCorners + 1];
	size_t numSideCorners = rRefElem.num(dim-1, side, 0);
	for (size_t co = 0; co < numSideCorners; ++co)
		vSideCorner[co] = vCornerCoords[rRefElem.id(dim-1, side, 0, co)];
	// we need another point if dim != TWorldDim:
	// take the highest-numbered corner of the next side, since
	// it is always different from the other points (is it not?)
	if constexpr (dim != TWorldDim)
	{
		vSideCorner[numSideCorners] =
			vCornerCoords[rRefElem.id(dim-1, (side+1)%rRefElem.num(dim-1), 0,
					  	  rRefElem.num(dim-1, (side+1)%rRefElem.num(dim-1), 0)-1)];
	}

//	Get the normal:
	ElementNormal<TWorldDim>(sideRoid, normalOut, vSideCorner);
//	Note: We assume that the for the standard ordering, the last line computes
//	the outer normal.
}

///	Computation of the side normal for a generic reference element:
template <int dim>
inline void SideNormal(ReferenceObjectID roid, MathVector<dim>& normalOut, int side, const MathVector<dim>* vCornerCoords);

template <>
inline void SideNormal<1>(ReferenceObjectID roid, MathVector<1>& normalOut, int side, const MathVector<1>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_EDGE: SideNormal<ReferenceEdge,1>(normalOut, side, vCornerCoords); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 1.");
	}
}

template <>
inline void SideNormal<2>(ReferenceObjectID roid, MathVector<2>& normalOut, int side, const MathVector<2>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_EDGE: SideNormal<ReferenceEdge,2>(normalOut, side, vCornerCoords); return;
		case ROID_TRIANGLE: SideNormal<ReferenceTriangle,2>(normalOut, side, vCornerCoords); return;
		case ROID_QUADRILATERAL: SideNormal<ReferenceQuadrilateral,2>(normalOut, side, vCornerCoords); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 2.");
	}
}

template <>
inline void SideNormal<3>(ReferenceObjectID roid, MathVector<3>& normalOut, int side, const MathVector<3>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_EDGE: SideNormal<ReferenceEdge,3>(normalOut, side, vCornerCoords); return;
		case ROID_TRIANGLE: SideNormal<ReferenceTriangle,3>(normalOut, side, vCornerCoords); return;
		case ROID_QUADRILATERAL: SideNormal<ReferenceQuadrilateral,3>(normalOut, side, vCornerCoords); return;
		case ROID_TETRAHEDRON: SideNormal<ReferenceTetrahedron,3>(normalOut, side, vCornerCoords); return;
		case ROID_PYRAMID: SideNormal<ReferencePyramid,3>(normalOut, side, vCornerCoords); return;
		case ROID_PRISM: SideNormal<ReferencePrism,3>(normalOut, side, vCornerCoords); return;
		case ROID_HEXAHEDRON: SideNormal<ReferenceHexahedron,3>(normalOut, side, vCornerCoords); return;
		case ROID_OCTAHEDRON: SideNormal<ReferenceOctahedron,3>(normalOut, side, vCornerCoords); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 3.");
	}
}

///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// ElementSideRayIntersection
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// wrapper class to distinguish reference dimension
template <typename TRefElem, int TWorldDim, int TRefDim = TRefElem::dim>
struct ElementSideRayIntersectionWrapper
{
	static bool apply(	size_t& sideOut,
						MathVector<TWorldDim>& GlobalIntersectionPointOut,
						MathVector<TRefElem::dim>& LocalIntersectionPoint,
						const MathVector<TWorldDim>& From, const MathVector<TWorldDim>& Direction,
						bool bPositiv, const MathVector<TWorldDim>* vCornerCoords)
	{
		UG_THROW("Not implemented.");
		return false;
	}
};

// specialization for 2d
template <typename TRefElem>
struct ElementSideRayIntersectionWrapper<TRefElem, 2, 2>
{
	static bool apply(	size_t& sideOut,
						MathVector<2>& GlobalIntersectionPointOut,
						MathVector<TRefElem::dim>& LocalIntersectionPoint,
						const MathVector<2>& From, const MathVector<2>& Direction,
						bool bPositiv, const MathVector<2>* vCornerCoords)
	{
		static const TRefElem& rRefElem = Provider<TRefElem>::get();

		// reference dimension
		const int dim = TRefElem::dim;

		// parameters
		number bc = 0., t = 0.;
		size_t p0 = 0, p1 = 0;

		// find side
		for(sideOut = 0; sideOut < rRefElem.num(dim-1); ++sideOut)
		{
			// get corners
			p0 = rRefElem.id(dim-1, sideOut, 0, 0);
			p1 = rRefElem.id(dim-1, sideOut, 0, 1);

			// if matched: break
			if(RayLineIntersection2d(	GlobalIntersectionPointOut, bc, t,
										vCornerCoords[p0], vCornerCoords[p1],
										From, Direction))
			{
				if ((bPositiv && t >= 0.0) || (!bPositiv && t <= 0.0))
					break;
			}
		}
		// if not found
		if(sideOut >= rRefElem.num(dim-1))
			UG_THROW("ElementSideRayIntersection: no cut side found.");

		// Compute local intersection
		VecScaleAdd(LocalIntersectionPoint, bc, rRefElem.corner(p1), 1.-bc, rRefElem.corner(p0));

		// true if found
		return true;
	}
};

// specialization for 3d
template <typename TRefElem>
struct ElementSideRayIntersectionWrapper<TRefElem, 3, 3>
{
	static bool apply(	size_t& sideOut,
						MathVector<3>& GlobalIntersectionPointOut,
						MathVector<TRefElem::dim>& LocalIntersectionPoint,
						const MathVector<3>& From, const MathVector<3>& Direction,
						bool bPositiv, const MathVector<3>* vCornerCoords)
	{
		static const TRefElem& rRefElem = Provider<TRefElem>::get();

		// reference dimension
		constexpr int dim = TRefElem::dim;

		// parameters
		number bc0 = 0., bc1 = 0., t = 0.;
		size_t p0 = 0, p1 = 0, p2 = 0;

		// find side
		for(sideOut = 0; sideOut < rRefElem.num(dim-1); ++sideOut)
		{
			// get corners
			p0 = rRefElem.id(dim-1, sideOut, 0, 0);
			p1 = rRefElem.id(dim-1, sideOut, 0, 1);
			p2 = rRefElem.id(dim-1, sideOut, 0, 2);

			// if match: break
			if(RayTriangleIntersection(	GlobalIntersectionPointOut, bc0, bc1, t,
										vCornerCoords[p0], vCornerCoords[p1], vCornerCoords[p2],
										From, Direction))
			{
				if ((bPositiv && t >= 0.0) || (!bPositiv && t <= 0.0))
					break;
			}

			// second triangle (only if 4 corners)
			if(rRefElem.num(dim-1, sideOut, 0) == 3) continue;

			// get corner number 4
			p1 = rRefElem.id(dim-1, sideOut, 0, 3);

			// if matched: break
			if(RayTriangleIntersection(	GlobalIntersectionPointOut, bc0, bc1, t,
										vCornerCoords[p0], vCornerCoords[p1], vCornerCoords[p2],
										From, Direction))
			{
				if ((bPositiv && t >= 0.0) || (!bPositiv && t <= 0.0))
					break;
			}
		}

		// if not found
		if(sideOut >= rRefElem.num(dim-1))
			UG_THROW("ElementSideRayIntersection: no cut side found.");

		// Compute local intersection
		VecScaleAdd(LocalIntersectionPoint,
		            (1.-bc0-bc1), rRefElem.corner(p0),
					bc0, rRefElem.corner(p1),
					bc1, rRefElem.corner(p2));

		// true if found
		return true;
	}
};


///////////////////////////////////////////////////////////////
/// ElementSideRayIntersection
/**
 * This function computes the side of element, that is intersected
 * by a Ray given in Parameter form by 'From + c * Direction'. The
 * intersection is choose to be at positive parameter if bPositiv == true,
 * else the intersection at negative parameter is choosen.
 * Global and local coordinates of the intersection points are returned
 * as well as the reference element number of the side.
 *
 * \param[in]	From							Point of Ray
 * \param[in]	Direction						Direction of Ray
 * \param[in]	bPositiv						Flag, whether to search in positiv of negative direction
 * \param[in]	vCornerCoords					Vector of corner coordinates
 * \param[out]	sideOut							side of intersection
 * \param[out]	GlobalIntersectionPointOut		Intersection Point (global)
 * \param[out]	LocalIntersectionPoint			Intersection Point (local)
 */
template <typename TRefElem, int TWorldDim>
bool ElementSideRayIntersection(	size_t& sideOut,
									MathVector<TWorldDim>& GlobalIntersectionPointOut,
									MathVector<TRefElem::dim>& LocalIntersectionPoint,
									const MathVector<TWorldDim>& From, const MathVector<TWorldDim>& Direction,
									bool bPositiv, const MathVector<TWorldDim>* vCornerCoords)
{
	UG_ASSERT(VecTwoNorm(Direction) > 0, "Direction must be non-zero vector.");
	return ElementSideRayIntersectionWrapper<TRefElem, TWorldDim>::
			apply(sideOut, GlobalIntersectionPointOut, LocalIntersectionPoint,
					From, Direction, bPositiv, vCornerCoords);
}


///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////
// ElementSideRayIntersection
///////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////

// wrapper class to distinguish reference dimension
template <int TDim, int TWorldDim>
struct SCVFofSCVRayIntersectionWrapper
{
	static bool apply(	size_t& sideOut, number& bc,
						MathVector<TWorldDim>& GlobalIntersectionPointOut,
						MathVector<TDim>& LocalIntersectionPoint,
						const MathVector<TWorldDim>& From, const MathVector<TWorldDim>& Direction,
						bool bPositiv, const MathVector<TWorldDim>* vCornerCoords)
	{
		UG_THROW("Not implemented.");
		return false;
	}
};

// specialization for 2d
template <>
struct SCVFofSCVRayIntersectionWrapper<2, 2>
{
	static bool apply(	size_t& sideOut, number& bc,
						MathVector<2>& GlobalIntersectionPointOut,
						MathVector<2>& LocalIntersectionPoint,
						const MathVector<2>& From, const MathVector<2>& Direction,
						bool bPositiv, const MathVector<2>* vCornerCoords)
	{
		static const ReferenceQuadrilateral& rRefElem = Provider<ReferenceQuadrilateral>::get();

		// reference dimension
		static constexpr int dim = 2;

		// parameters
		bc = 0.;
		number t = 0.;
		size_t p0 = 0, p1 = 0;

		// find side
		for(sideOut = 0; sideOut < rRefElem.num(0); ++sideOut)
		{
			// get corners
			p0 = rRefElem.id(dim-1, sideOut, 0, 0);
			p1 = rRefElem.id(dim-1, sideOut, 0, 1);

			// if matched: break
			if(RayLineIntersection2d(	GlobalIntersectionPointOut, bc, t,
										vCornerCoords[p0], vCornerCoords[p1],
										From, Direction))
			{
			//	skip if one root-point scvf
				if(fabs(t) <= std::numeric_limits<number>::epsilon() * 10)
					continue;

			//	upwind / downwind switch
				if ((bPositiv && t >= 0.0) || (!bPositiv && t <= 0.0))
					break;
			}
		}

		// if not found
		if(sideOut >= rRefElem.num(0))
			UG_THROW("Side not found.");

		// Compute local intersection
		VecScaleAdd(LocalIntersectionPoint, bc, rRefElem.corner(p1), 1.-bc, rRefElem.corner(p0));

		//	parameter of intersection should loop always from center to edgeMidpoint
		//	thus, for side 2 this is correct, but for side 1 we have to invert direction
		if(sideOut == 1) bc = 1. - bc;

		// true if found on scvf
		if(sideOut == 1 || sideOut == 2) return true;
		// false if found on element corner
		else return false;
	}
};

///////////////////////////////////////////////////////////////
// SCVFofSCVRayIntersection
/**
 * This function computes if another scvf of a scv is intersected
 * by a Ray given in Parameter form by 'Root + c * Direction'. The
 * intersection is choose to be at positive parameter if bPositiv == true,
 * else the intersection at negative parameter is choosen.
 * Global and local coordinates of the intersection points are returned
 * as well as the local number of the side of the scv that matches the intersection.
 *
 * \param[in]	Root							Point of Ray
 * \param[in]	Direction						Direction of Ray
 * \param[in]	bPositiv						Flag, whether to search in positiv of negative direction
 * \param[in]	vCornerCoords					Vector of corner coordinates
 * \param[out]	sideOut							side of intersection
 * \param[out]	bc								line parameter in [0,1] indicating position
 * 												of intersection point on scvf in direction
 * 												center to edge
 * \param[out]	GlobalIntersectionPointOut		Intersection Point (global)
 * \param[out]	LocalIntersectionPoint			Intersection Point (local)
 * \returns true if intersected with another scvf of the scv, false else
 */
template <int TDim, int TWorldDim>
bool SCVFofSCVRayIntersection(	size_t& sideOut, number& bc,
                              	MathVector<TWorldDim>& GlobalIntersectionPointOut,
                              	MathVector<TDim>& LocalIntersectionPoint,
                              	const MathVector<TWorldDim>& Root, const MathVector<TWorldDim>& Direction,
                              	bool bPositiv, const MathVector<TWorldDim>* vCornerCoords)
{
	UG_ASSERT(VecTwoNorm(Direction) > 0, "Direction must be non-zero vector.");
	return SCVFofSCVRayIntersectionWrapper<TDim, TWorldDim>::
			apply(sideOut, bc, GlobalIntersectionPointOut, LocalIntersectionPoint,
					Root, Direction, bPositiv, vCornerCoords);
}



///////////////////////////////////////////////////////////////
/// Extension of an element (in a given dimension)
///////////////////////////////////////////////////////////////


template <typename TVector>
inline void ComputeElementExtensionsSqForEdges(const TVector* vCornerCoords, TVector &ext)
{
	VecElemProd(ext, vCornerCoords[0], vCornerCoords[1]);
}


template <int TWorldDim, int ncorners>
inline void ComputeElementExtensionsSq(const MathVector<TWorldDim>* vCornerCoords, MathVector<TWorldDim> &ext)
{
	// compute center
	MathVector<TWorldDim> mid = 0.0;
	for (int i=ncorners-1; i>=0; --i){
		VecAppend(mid, vCornerCoords[i]);
	}
	VecScale(mid, mid, 1.0/ncorners);

	// compute
	ext = 0.0;
	MathVector<TWorldDim> aux;
	for (int i=ncorners-1; i>=0; --i){
		VecSubtract(aux, mid, vCornerCoords[i]);
		VecElemProd(aux, aux, aux);
		VecAppend(ext, aux);
	}
	VecScale(ext, ext, 1.0/ncorners);
}

///////////////////////////////////////////////////////////////
//	extensions of an element
///////////////////////////////////////////////////////////////


template <int dim>
inline void ElementExtensionsSq(ReferenceObjectID roid, MathVector<dim>& ext, const MathVector<dim>* vCornerCoords);

template <>
inline void ElementExtensionsSq<1>(ReferenceObjectID roid, MathVector<1>& ext, const MathVector<1>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: 			ext = 0; return;
		case ROID_EDGE: 			ComputeElementExtensionsSqForEdges(vCornerCoords,ext); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 1.");
	}
}

template <>
inline void ElementExtensionsSq<2>(ReferenceObjectID roid, MathVector<2>& ext, const MathVector<2>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: 			ext = 0; return;
		case ROID_EDGE: 			ComputeElementExtensionsSqForEdges(vCornerCoords,ext); return;
		case ROID_TRIANGLE: 		ComputeElementExtensionsSq<2,3>(vCornerCoords,ext); return;
		case ROID_QUADRILATERAL: 	ComputeElementExtensionsSq<2,4>(vCornerCoords,ext); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 2.");
	}
}

template <>
inline void ElementExtensionsSq<3>(ReferenceObjectID roid, MathVector<3>& ext, const MathVector<3>* vCornerCoords)
{
	switch(roid)
	{
		case ROID_VERTEX: 			ext = 0; return;
		case ROID_EDGE: 			ComputeElementExtensionsSqForEdges(vCornerCoords,ext); return; 
		case ROID_TRIANGLE: 		ComputeElementExtensionsSq<3,3>(vCornerCoords,ext); return;
		case ROID_QUADRILATERAL: 	ComputeElementExtensionsSq<3,4>(vCornerCoords,ext); return;
		case ROID_TETRAHEDRON: 		ComputeElementExtensionsSq<3,4>(vCornerCoords,ext); return;
		case ROID_PYRAMID: 			ComputeElementExtensionsSq<3,5>(vCornerCoords,ext); return;
		case ROID_PRISM: 			ComputeElementExtensionsSq<3,6>(vCornerCoords,ext); return;
		case ROID_HEXAHEDRON: 		ComputeElementExtensionsSq<3,8>(vCornerCoords,ext); return;
		case ROID_OCTAHEDRON: 		ComputeElementExtensionsSq<3,6>(vCornerCoords,ext); return;
		default: UG_THROW("ReferenceObject "<<roid<<" not found in dim 3.");
	}
}


} // end namespace ug

#endif
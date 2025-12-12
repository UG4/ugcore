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

#ifndef __H__UGMATH__MATH_UTIL__
#define __H__UGMATH__MATH_UTIL__

//#include <cstdlib>
#include "math_constants.h"
#include "../math_vector_matrix/math_vector.h"
#include "../math_vector_matrix/math_vector_functions.h"
#include "../ugmath_types.h"
#include "../../ug_config.h"
//#include "eigenvalues.h"

// this is quite hacky...
// some environments (on windows) seem to define a small type or variable.
// this causes errors with parameter names etc. It is thus undefined below.
#ifdef small
	#undef small
#endif

namespace ug {

/**
 * \defgroup ugbase_mathutil Math Utilities
 * \ingroup ugbase_math
 * \{
 */

////////////////////////////////////////////////////////////////////////
//	deg_to_rad
template <typename TNumber>
inline TNumber
deg_to_rad(TNumber deg);

////////////////////////////////////////////////////////////////////////
//	rad_to_deg
template <typename TNumber>
inline TNumber
rad_to_deg(TNumber rad);

////////////////////////////////////////////////////////////////////////
//	urand
///	uniform distributed random numbers in [lowerBound, upperBound[. Use srand to set a seed.
template <typename TNumber>
TNumber
urand(TNumber lowerBound, TNumber upperBound);

////////////////////////////////////////////////////////////////////////
//	clip
///	clips a number to the given interval [lowerBound, upperBound].
template <typename TNumber>
TNumber
clip(TNumber val, TNumber lowerBound, TNumber upperBound);

////////////////////////////////////////////////////////////////////////
///	returns the square of a value (val*val)
template <typename TNumber>
inline TNumber sq(TNumber val);

////////////////////////////////////////////////////////////////////////
///	finds a normal to the given vector in 3d.
UG_API
bool FindNormal(vector3& normOut, const vector3& v);

////////////////////////////////////////////////////////////////////////
///	constructs a orthonormal matrix given a vector in 3d.
/**
 * Given a vector in 3d, this method constructs an orthonormal system.
 * The given vector correspond to one of the matrix columns, the others
 * are found so that matOut is orthonormal (inv(matOut) == trans(matOut)).
 * Please note that matOut is not unique.
 *
 * This method assumes that the matrix is accessed as follows:
 * m[colInd][rowInd].
 *
 * \param matOut: Matrix
 * \param v: vector
 * \param vColInd: specifies in which column to put v.
 */
UG_API
bool ConstructOrthonormalSystem(matrix33& matOut, const vector3& v,
								size_t vColInd);

////////////////////////////////////////////////////////////////////////
///	calculates the center of a point-set
template <typename vector_t>
void CalculateCenter(vector_t& centerOut, const vector_t* pointSet,
					 size_t numPoints);

////////////////////////////////////////////////////////////////////////
///	Calculates the barycenter of a triangle (1/3) * (p1+p2+p3)
template <typename vector_t>
vector_t
TriangleBarycenter(const vector_t& p1, const vector_t& p2, const vector_t& p3);

////////////////////////////////////////////////////////////////////////
///	Calculates the circumcenter of a triangle
/**	Calculates the center of the circle, on which all three points lie.
 * In cases where this is not possible (all three points lie on a line), the
 * method returns false.
 * If the method succeeds, the calculated center is written to centerOut.
 * \{ */
UG_API
bool TriangleCircumcenter(vector2& centerOut, const vector2& p1,
						  const vector2& p2, const vector2& p3);

UG_API
bool TriangleCircumcenter(vector3& centerOut, const vector3& p1,
						  const vector3& p2, const vector3& p3);
/**	\} */

////////////////////////////////////////////////////////////////////////
///	Calculates the covariance matrix of a given point-set.
/**	
 * Please note that you have to specify the point-set together with
 * its center-point.
 */
UG_API
void CalculateCovarianceMatrix(matrix33& matOut, const vector3* pointSet,
							  const vector3& center, size_t numPoints);

////////////////////////////////////////////////////////////////////////
///	Finds the plane which minimizes the distance to the given point set.
/**	The method returns the true if the calculation was a success and false
 * if not. If true is returned, then centerOut and normalOut will contain
 * the center and the normal of the plane which minimizes the distance to
 * the given point set.
 */
UG_API
bool FindClosestPlane(vector3& centerOut, vector3& normalOut,
					  const vector3* pointSet, size_t numPoints);

////////////////////////////////////////////////////////////////////////
///	transforms points from 3d space to 2d space.
/**
 * This method calculates the plane that minimizes the distances
 * between the points and the plane. It then transforms all points
 * to 2d so that the normal of the new min-distance-plane would
 * point along the z-axis. It the projects all points onto
 * the x-y-plane using the trivial projection.
 *
 * \param pointSetOut: An array of 2d-points of size numPoints.
 * \param pointSet: An array of 3d-points of size numPoints.
 * \param numPoints: Specifies the size of the point sets.
 */
UG_API
bool TransformPointSetTo2D(vector2* pointSetOut, const vector3* pointSet,
						  size_t numPoints);

////////////////////////////////////////////////////////////////////////
///	finds the projection of v onto the line defined by v0 and v1
/**
 * projects v onto the line defined by v0 and v1.
 * The projected point is returned in vOut.
 *
 * returns s so that vOut = (1.0-s)*v0 + s*v1
 */
template <typename vector_t>
number DropAPerpendicular(vector_t& vOut, const vector_t& v,
							const vector_t& v0, const vector_t& v1);


////////////////////////////////////////////////////////////////////////
///	returns the point described by a relative ray coordinate
template <typename vector_t>
vector_t PointOnRay(const vector_t& from, const vector_t& dir, number s);

////////////////////////////////////////////////////////////////////////
///	finds the projection of v onto the ray defined by from and dir
/**
 * projects v onto the ray defined by from and dir.
 * The projected point is returned in vOut.
 *
 * returns s so that vOut = from + s*dir
 */
template <typename vector_t>
number ProjectPointToRay(vector_t& vOut, const vector_t& v,
							const vector_t& from, const vector_t& dir);


////////////////////////////////////////////////////////////////////////
///	finds the projection of v onto the line defined by from and to
/** projects v onto the line defined by from and to.
 * The projected point is returned in vOut.
 *
 * returns s so that vOut = from + s * (to-from)*/
template <typename vector_t>
number ProjectPointToLine(vector_t& vOut, const vector_t& v,
						  const vector_t& from, const vector_t& to);


////////////////////////////////////////////////////////////////////////
///	calculates the distance of a point to a line segment
/*
 * the coordinates of the point are given through v.
 * the line-segment is defined by its endpoints v1 and v2.
 *
 * The method returns the distance.
 */
template <typename vector_t>
inline
number DistancePointToLine(const vector_t& v, const vector_t& v1,
						  const vector_t& v2);

template <typename vector_t>
number DistancePointToLine(number& tOut, const vector_t& v,
						   const vector_t& v1, const vector_t& v2);
						  
////////////////////////////////////////////////////////////////////////
///	calculates the distance of a point to a ray
/*
 * the coordinates of the point are given through v.
 * the ray is defined by from (an arbitrary point on the ray) and
 * dir (the direction of the ray).
 *
 * The method returns the distance.
 */
template <typename vector_t>
inline
number DistancePointToRay(const vector_t& v, const vector_t& from,
						  const vector_t& dir);

template <typename vector_t>
inline
number DistancePointToRay(number& tOut, const vector_t& v,
						  const vector_t& from, const vector_t& dir);

////////////////////////////////////////////////////////////////////////
///	calculates the minimal distance of a point to a triangle.				
/**	Distance is measured between the point p and the point on the
 *	triangle which is closest to p.
 *
 *	The triangle is defined by the three corner-points v1, v2 and v3.
 *
 *	n specifies the normal of the triangle. It does not have do be normalized.
 *
 *	The closest point on the triangle is written to vOut. Its barycentric
 *	coordinate is written to bc1Out and bc2Out.
 *
 *	The method returns the minimal distance of the point to the triangle.
 */
template <typename vector_t>
number DistancePointToTriangle(vector_t& vOut, number& bc1Out, number& bc2Out,
							const vector_t& p, const vector_t& v1, const vector_t& v2,
							const vector_t& v3, const vector_t& n);

////////////////////////////////////////////////////////////////////////
///	Calculates the distance between the specified point and the plane.
template <typename vector_t>
number DistancePointToPlane(const vector_t& v, const vector_t& p,
							const vector_t& n);

////////////////////////////////////////////////////////////////////////
///	projects v onto the plane defined by the point p and the planes normal n.
/**	The result is written to vOut.*/
template <typename vector_t>
void ProjectPointToPlane(vector_t& vOut, const vector_t& v,
						const vector_t& p, const vector_t& n);

////////////////////////////////////////////////////////////////////////
//	RayPlaneIntersection
///	calculates the intersection of the ray rayFrom+t*rayDir and the plane (x-p)*n=0.
template <typename vector_t>
bool RayPlaneIntersection(vector_t& vOut, number& tOut,
						  const vector_t& rayFrom, const vector_t& rayDir,
						  const vector_t& p, const vector_t& n);

////////////////////////////////////////////////////////////////////////
//	RayRayIntersection2d
///	calculates the intersection of two Rays in 2d
/** If the two rays intersect (are not parallel), the method returns true.
 * It writes its results to the parameters vOut (intersection point),
 * t0Out (vOut = p0 + t0Out * dir0) and t1Out (vOut = p1 + t1Out * dir1).
 *
 * You have to specify a point of each ray (p0 and p1) and the directions of
 * each ray (dir0 and dir1).
 */
template <typename vector_t>
bool RayRayIntersection2d(vector_t &vOut, number& t0Out, number& t1Out,
						   const vector_t &p0, const vector_t &dir0,
						   const vector_t &p1, const vector_t &dir1);

////////////////////////////////////////////////////////////////////////
//	RayLineIntersection2d
///	calculates the intersection of a ray with a Line in 2d
/**
 * You have to pass the line corners through p0 and p1
 * together with a point on the ray (vFrom) and the rays
 * direction (vDir).
 *
 * If the method succeeds (the ray intersects the line)
 * the methods returns true and writes the position of
 * the intersection to vOut.
 * Furthermore the local (barycentric) coordinate of the
 * intersection is written to bcOut. tOut
 * will contain the local coordinate of the intersection
 * regarding the rays parameter form.
 */
template <typename vector_t>
bool RayLineIntersection2d(vector_t &vOut, number& bcOut, number& tOut,
						   const vector_t &p0, const vector_t &p1,
						   const vector_t &vFrom, const vector_t &vDir,
						   number sml = 0);


////////////////////////////////////////////////////////////////////////
//	LineLineIntersection2d
///	calculates the intersection of a two lines in 2d
/**
 * The lines are specified through their endpoints:
 * (from0, to0) describes the first line, (from1, to1) the second line.
 *
 * If the method succeeds (the lines intersects)
 * the methods returns true and writes the position of
 * the intersection to vOut.
 * t0Out and t1Out will contain the local coordinate of the intersection
 * regarding the lines parameter form.
 */
template <typename vector_t>
bool LineLineIntersection2d(vector_t &vOut, number& t0Out, number& t1Out,
						   const vector_t &from0, const vector_t &to0,
						   const vector_t &from1, const vector_t &to1,
						   const number threshold = 0);


////////////////////////////////////////////////////////////////////////
///	intersects two 3d rays
/**	Returns true, if the rays intersect. If they don't, the output
 * parameters are filled with the points, at which the rays have the
 * minimal distance.
 * Specify each ray by one point on the ray (p0 and p1) and the direction
 * of the ray (dir0 and dir1).
 *
 * Output-parameters: aOut, bOut represent the closest points on the two rays.
 *
 * This method internally uses the IntersectLineSegments algorithm by Graham Rhodes.
 * Please have a look at lineintersect_utils.h for more information.*/
UG_API
bool RayRayIntersection3d(vector3& aOut, vector3& bOut,
						  const vector3& p0, const vector3& dir0,
						  const vector3& p1, const vector3& dir1);

////////////////////////////////////////////////////////////////////////
///	intersects two 3d line segments
/**	Returns true, if the lines intersect. If they don't, the output
 * parameters are filled with the points, at which the lines have the
 * minimal distance.
 * Specify lines a (given by endpoints a1, a2) and b (given by endpoints b1, b2).
 *
 * Output-parameters: aOut, bOut represent the closest points on the two lines.
 *
 * This method internally uses the IntersectLineSegments algorithm by Graham Rhodes.
 * Please have a look at lineintersect_utils.h for more information.*/
UG_API
bool LineLineIntersection3d(vector3& aOut, vector3& bOut,
							const vector3& a1, const vector3& a2,
						  	const vector3& b1, const vector3& b2);

////////////////////////////////////////////////////////////////////////
///	Calculates the parameter values at wich two rays are closest.
/**	Returns true if the calculation was successfull (e.g. if the rays were not parallel).
 * \sa RayRayIntersection2d, RayRayIntersection3d, LineLineProjection*/
template <typename vector_t>
bool RayRayProjection(number& t1Out, number& t2Out,
						const vector_t& from1, const vector_t& dir1,
						const vector_t& from2, const vector_t& dir2);

////////////////////////////////////////////////////////////////////////
///	calculates the closest point between the rays through the given lines.
/**	returns true if the closest points lie inside the line segments.*/
template <typename vector_t>
bool LineLineProjection(number& t1Out, number& t2Out,
						  const vector_t& a1, const vector_t& a2,
						  const vector_t& b1, const vector_t& b2);

////////////////////////////////////////////////////////////////////////
///	calculates the distance between two 3d line segments
/**	Calculates the distance between the two finite line segments
 * a (given by endpoints a1, a2) and b (given by endpoints b1, b2).
 *
 * This method internally uses the IntersectLineSegments algorithm by Graham Rhodes.
 * Please have a look at lineintersect_utils.h for more information.*/
UG_API
number DistanceLineToLine(const vector3& a1, const vector3& a2,
						  const vector3& b1, const vector3& b2);

////////////////////////////////////////////////////////////////////////
//	RayTriangleIntersection
///	calculates the intersection of a ray with a triangle
/**
 * vector_t has to feature a x, y and z component.
 *
 * You have to pass the triangles corners through p0, p1, p2
 * together with a point on the ray (vFrom) and the rays
 * direction (vDir).
 *
 * If the method succeeds (the ray intersects the triangle)
 * the methods returns true and writes the position of
 * the intersection to vOut.
 * Furthermore the local (barycentric) coordinates of the
 * intersection are written to bc1Out and bc2Out. tOut
 * will contain the local coordinate of the intersection
 * regarding the rays parameter form.
 */
template <typename vector_t>
bool RayTriangleIntersection(vector_t &vOut, number& bc1Out, number& bc2Out, number& tOut,
						   const vector_t &p0, const vector_t &p1, const vector_t &p2, 
						   const vector_t &vFrom, const vector_t &vDir,
						   const number small = SMALL);

////////////////////////////////////////////////////////////////////////
//	RayTriangleIntersection
///	calculates the intersection of a ray with a triangle
/**
 * vector_t has to feature a x, y and z component.
 *
 * You have to pass the triangles corners through p0, p1, p2
 * together with a point on the ray (vFrom) and the rays
 * direction (vDir).
 *
 * If the method succeeds (the ray intersects the triangle)
 * the methods returns true and writes the position of
 * the intersection to vOut.
 */
template <typename vector_t> inline
bool RayTriangleIntersection(vector_t &vOut, const vector_t &p0,
						   const vector_t &p1, const vector_t &p2, 
						   const vector_t &vFrom, const vector_t &vDir);

////////////////////////////////////////////////////////////////////////
//	RayBoxIntersection
///	checks if a ray is intersecting a box
/**
 * vector_t has to feature a static 'Size' member and index access
 *
 * \param rayFrom:
 * \param rayDir:
 * \param boxMin:
 * \param boxMax:
 * \param tNearOut: can be nullptr
 * \param tFarOut: can be nullptr
 *
 * \sa LineBoxIntersection
 */
template <typename vector_t>
bool RayBoxIntersection(const vector_t& rayFrom, const vector_t& rayDir,
						const vector_t& boxMin, const vector_t& boxMax,
						number* tNearOut = nullptr, number* tFarOut = nullptr);

////////////////////////////////////////////////////////////////////////
///	checks whether the given line-segment (v1, v2) intersect the given box.
/**
 * vector_t has to feature a static 'Size' member and index access
 *
 * \sa RayBoxIntersection
 */
template <typename vector_t>
bool LineBoxIntersection(const vector_t& v1, const vector_t& v2,
						const vector_t& boxMin, const vector_t& boxMax);

////////////////////////////////////////////////////////////////////////
///	checks whether the line segment (v1, v2) intersect the given sphere.
/**	The function returns the number of intersections and the values of s at
 * which the line segment intersects.
 */
template <typename vector_t>
int RaySphereIntersection(number& s1Out, number& s2Out,
						  const vector_t& v, const vector_t& dir,
						  const vector_t& center, number radius);

template <typename vector_t>
int LineSphereIntersection(number& s1Out, number& s2Out,
						  const vector_t& v1, const vector_t& v2,
						  const vector_t& center, number radius);


////////////////////////////////////////////////////////////////////////
///	returns the parameter values at which a given ray intersects an infinite cylinder.
UG_API
bool RayCylinderIntersection(number& tMinOut, number& tMaxOut, const vector3& rayFrom,
							 const vector3& rayDir, const vector3& cylCenter,
							 const vector3& cylAxis, number cylRadius);

////////////////////////////////////////////////////////////////////////
//	TriangleTriangleIntersection
///	checks whether two triangles intersect and returns the intervals, if they do.
/**	Uses code by Tomas Moller.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * returns true if the triangles intersect. If ip1Out and ip2Out are specified
 * and if the triangles intersect, ip1Out and ip2Pit will contain the endpoints
 * of the line segment which resembles the intersection. Please specify either both
 * or none.
 */
UG_API
bool TriangleTriangleIntersection(const MathVector<3>& p0, const MathVector<3>& p1,
								  const MathVector<3>& p2, const MathVector<3>& q0,
								  const MathVector<3>& q1, const MathVector<3>& q2,
								  MathVector<3>* ip1Out = nullptr,
								  MathVector<3>* ip2Out = nullptr,
								  number snapThreshold = SMALL);

////////////////////////////////////////////////////////////////////////
//	TriBoxIntersection
//	opposed to most of the other methods here, TriBoxIntersection is not
//	a template method. Instead it directly operates on the built in
//	vector3 type. This is due to its implementation, which is taken
//	from 'Game Programming Gems'. Please have a look at the copyright notice
//	in 'tri_box.cpp'.
////////////////////////////////////////////////////////////////////////
///	checks whether a triangle and an axis-aligned box intersect.
/**
 * \param p0, p1, p2: the corners of the triangle.
 * \param boxMin: has to contain the minimal-coordinates of the box.
 * \param boxMax: has to contain the maximal-coordinates of the box.
 */
UG_API
bool TriangleBoxIntersection(const MathVector<3>& p0, const MathVector<3>& p1,
							 const MathVector<3>& p2,
							 const MathVector<3>& boxMin, const MathVector<3>& boxMax);
						
////////////////////////////////////////////////////////////////////////
///	checks whether two boxes intersect.
/**	Make sure that all parameters have the same size.*/
template <typename vector_t>
bool BoxBoxIntersection(const vector_t& box1Min, const vector_t& box1Max,
						const vector_t& box2Min, const vector_t& box2Max);

////////////////////////////////////////////////////////////////////////
//	BoxBoundProbe
///	Returns true if the point lies inside or on the boundary of the box
/**	v.size() has to be the same or less than boxMin.size() and boxMax.size().
 *	If this condition is met, the method supports vectors of arbitrary dimension.*/
template <typename vector_t>
bool BoxBoundProbe(const vector_t& v, const vector_t& boxMin,
					const vector_t& boxMax);

////////////////////////////////////////////////////////////////////////
///	calculates the are of the triangle defined by p1, p2 and p3
template <typename vector_t>
number TriangleArea(const vector_t& p1, const vector_t& p2,
					const vector_t& p3);

////////////////////////////////////////////////////////////////////////
///	calculates the are of the triangle defined by p1, p2 and p3
template <typename vector_t>
number QuadrilateralArea(const vector_t& p1, const vector_t& p2,
						 const vector_t& p3, const vector_t& p4);

////////////////////////////////////////////////////////////////////////
///	the returned degree lies between 0 and 1. The closer to 1 the better.
/**
 * returns the minimal dot-product of each normal of the triangle-corners
 * with the triangles normal.
 */
template <typename vector_t>
number GeometricApproximationDegree(vector_t& n1, vector_t& n2, vector_t& n3,
									vector_t& tn);

////////////////////////////////////////////////////////////////////////
///	returns a value between 0 and 1. The higher the better.
/**	The quality is measured by comparing the triangles area with the
 *	sum of the squared length of each side.
 *
 *	a triangle whose sides have all zero length is considered to be
 *	a bad triangle - 0 is returned.*/
template <typename vector_t>
number TriangleQuality_Area(const vector_t& p1, const vector_t& p2,
							const vector_t& p3);

////////////////////////////////////////////////////////////////////////
///	Returns true if the point lies inside or on the boundary of a triangle
/**
 * This method should only be used if v, v0, v1, and v2 do lie in the same
 * x-y-plane.
 */
template <typename vector_t>
bool PointIsInsideTriangle(const vector_t& v, const vector_t& v0,
							const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////////////
///	Returns true if the point lies inside or on a side of the given triangle
/**
 * This method should only be used if v, v0, v1, and v2 do lie in the same
 * x-y-plane.
 *
 * The method has a higher accuracy than PointIsInsideTriangle, however it also
 * has a higher performance overhead.
 */
template <typename vector_t>
bool PointIsInsideTriangle_HighAcc(const vector_t& v, const vector_t& v0,
								   const vector_t& v1, const vector_t& v2);

////////////////////////////////////////////////////////////////////////
///	Returns true if the point lies inside or on the boundary of a quadrilateral
/**
 * This method should only be used if v, v0, v1, v2, and v3 do lie in the same
 * x-y-plane.
 */
template <typename vector_t>
bool PointIsInsideQuadrilateral(const vector_t& v, const vector_t& v0,
								const vector_t& v1, const vector_t& v2,
								const vector_t& v3);

////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron
///	Returns true if the point lies inside or on the boundary of a tetrahedron
/**
 * This method does not care about the orientation of the tetrahedron.
 * This however makes it a little slower than a method that only
 * works for correctly orientated tetrahedrons.
 */
template <typename vector_t>
bool PointIsInsideTetrahedron(const vector_t& v, const vector_t& v0, const vector_t& v1,
							  const vector_t& v2, const vector_t& v3);

////////////////////////////////////////////////////////////////////////
//	CalculateTetrahedronVolume - mstepnie
UG_API
number CalculateTetrahedronVolume(const vector3& a, const vector3& b,
								  const vector3& c, const vector3& d);

//	CalculatePyramidVolume - mscherer
number CalculatePyramidVolume(const vector3& v0, const vector3& v1,
		const vector3& v2, const vector3& v3, const vector3& v4);

//	CalculatePrismVolume - mscherer
number CalculatePrismVolume(const vector3& v0, const vector3& v1,
		const vector3& v2, const vector3& v3, const vector3& v4,
		const vector3& v5);

//	CalculateHexahedronVolume - mscherer
number CalculateHexahedronVolume(const vector3& v0, const vector3& v1,
		const vector3& v2, const vector3& v3, const vector3& v4,
		const vector3& v5, const vector3& v6, const vector3& v8);

//	CalculateHexahedronVolume - mstepnie
/**
 * (compare to lib_disc/common/geometry_util.h)
 * This function returns the volume of an octhedron in 3d
 * by calculating the volumes of the upper and lower pyramid
 * the octahedron consists of.
 * The pyramidal volumes are computed via: V = 1/3 * (S * h) with
 * - S is the area of the base
 * - h is the height
 *
 * The corner coordinates must be given as prescribed by the reference element
 *
 * \param[in]	a,...,f			3d Vectors of corner coordinates (6 corners)
 * \return 		number			Volume of Octahedron
 */
number CalculateOctahedronVolume(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e,
		const vector3& f);

////////////////////////////////////////////////////////////////////////
//	BinomialCoefficient
///	Returns the BinomialCoefficient
/**
 * Returns:
 * 	  n!
 * ---------
 * k! (n-k)!
 */
UG_API
int BinomCoeff(int n, int k);

////////////////////////////////////////////////////////////////////////
//	ReflectVectorAtPlane
///	reflects a vector at a plane
/**
 * This method reflects a given vector at a plane specified by n and r0, i.e.
 * the plane P is given such that for all \f$ \vec{x} in P \f$ it holds that
 * \f[
 * 	(\vec{x} - \vec{r0}) \cdot \vec{n} = 0.
 * \f]
 *
 * @param vReflectedOut		the reflected vector
 * @param v					the original vector, that is reflected
 * @param n					normal for plane specification
 * @param r0				root point of plane
 */
template <typename vector_t>
void ReflectVectorAtPlane(vector_t& vReflectedOut, const vector_t& v,
                          const vector_t& n, const vector_t& r0);

// end group ugbase_math
/// \}

}//	end of namespace

////////////////////////////////////////////////
//	include implementation
#include "math_util_impl.hpp"

#endif

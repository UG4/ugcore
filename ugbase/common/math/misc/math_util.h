//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d14

#ifndef __H__UGMATH__MATH_UTIL__
#define __H__UGMATH__MATH_UTIL__

#include <cstdlib>
#include "../math_vector_matrix/math_vector.h"
#include "../math_vector_matrix/math_vector_functions.h"
#include "eigenvalues.h"

namespace ug
{

const number SMALL = 1.0e-12;
const number PI = 3.14159265;

////////////////////////////////////////////////////////////////////////
//	deg_to_rad
template <class TNumber>
inline TNumber
deg_to_rad(TNumber deg);

////////////////////////////////////////////////////////////////////////
//	rad_to_deg
template <class TNumber>
inline TNumber
rad_to_deg(TNumber rad);

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
///	finds a normal to the given vector in 3d.
bool FindNormal(ug::vector3& normOut, const ug::vector3& v);

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
bool ConstructOrthonormalSystem(ug::matrix33& matOut, const ug::vector3& v,
								size_t vColInd);

////////////////////////////////////////////////////////////////////////
///	calculates the center of a point-set
template <class vector_t>
void CalculateCenter(vector_t& centerOut, const vector_t* pointSet,
					 size_t numPoints);

////////////////////////////////////////////////////////////////////////
///	Calculates the covariance matrix of a given point-set.
/**	
 * Please note that you have to specify the point-set together with
 * its center-point.
 */
void CalculateCovarianceMatrix(matrix33& matOut, const ug::vector3* pointSet,
							  const ug::vector3& center, size_t numPoints);

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
bool TransformPointSetTo2D(ug::vector2* pointSetOut, const ug::vector3* pointSet,
						  size_t numPoints);

////////////////////////////////////////////////////////////////////////
///	finds the projection of v onto the line defined by v0 and v1
/**
 * projects v onto the line defined by v0 and v1.
 * The projected point is returned in vOut.
 *
 * returns s so that vOut = (1.0-s)*v0 + s*v1
 */
template <class vector_t>
number DropAPerpendicular(vector_t& vOut, const vector_t& v,
							const vector_t& v0, const vector_t& v1);

////////////////////////////////////////////////////////////////////////
///	finds the projection of v onto the ray defined by from and dir
/**
 * projects v onto the ray defined by from and dir.
 * The projected point is returned in vOut.
 *
 * returns s so that vOut = from + s*dir
 */
template <class vector_t>
number ProjectPointToRay(vector_t& vOut, const vector_t& v,
							const vector_t& from, const vector_t& dir);

////////////////////////////////////////////////////////////////////////
///	calculates the distance of a point to a line segment
/*
 * the coordinates of the point are given through v.
 * the line-segment is defined by its endpoints v1 and v2.
 *
 * The method returns the distance.
 */
template <class vector_t>
inline
number DistancePointToLine(const vector_t& v, const vector_t& v1,
						  const vector_t& v2);

template <class vector_t>
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
template <class vector_t>
inline
number DistancePointToRay(const vector_t& v, const vector_t& from,
						  const vector_t& dir);

template <class vector_t>
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
template <class vector_t>
number DistancePointToTriangle(vector_t& vOut, number& bc1Out, number& bc2Out,
							const vector_t& p, const vector_t& v1, const vector_t& v2,
							const vector_t& v3, const vector_t& n);

////////////////////////////////////////////////////////////////////////
///	projects v onto the plane defined by the point p and the planes normal n.
/**	The result is written to vOut.*/
template <class vector_t>
void ProjectPointToPlane(vector_t& vOut, const vector_t& v,
						const vector_t& p, const vector_t& n);


////////////////////////////////////////////////////////////////////////
//	RayPlaneIntersection
///	calculates the intersection of the ray rayFrom+t*rayDir and the plane (x-p)*n=0.
template <class vector_t>
bool RayPlaneIntersection(vector_t& vOut, number& tOut,
						  const vector_t& rayFrom, const vector_t& rayDir,
						  const vector_t& p, const vector_t& n);

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
template <class vector_t>
bool RayLineIntersection2d(vector_t &vOut, number& bcOut, number& tOut,
						   const vector_t &p0, const vector_t &p1,
						   const vector_t &vFrom, const vector_t &vDir);

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
template <class vector_t>
bool RayTriangleIntersection(vector_t &vOut, number& bc1Out, number& bc2Out, number& tOut,
						   const vector_t &p0, const vector_t &p1, const vector_t &p2, 
						   const vector_t &vFrom, const vector_t &vDir);

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
template <class vector_t> inline
bool RayTriangleIntersection(vector_t &vOut, const vector_t &p0,
						   const vector_t &p1, const vector_t &p2, 
						   const vector_t &vFrom, const vector_t &vDir);

////////////////////////////////////////////////////////////////////////
//	RayBoxIntersection
///	checks if a ray is intersecting a box
/**
 * vector_t has to feature a x, y and z component.
 *
 * \param rayFrom:
 * \param rayDir:
 * \param boxMin:
 * \param boxMax:
 * \param tNearOut: can be NULL
 * \param tFarOut: can be NULL
 *
 * \sa LineBoxIntersection
 */
template <class vector_t>
bool RayBoxIntersection(const vector_t& rayFrom, const vector_t& rayDir,
						const vector_t& boxMin, const vector_t& boxMax,
						number* tNearOut = NULL, number* tFarOut = NULL);

////////////////////////////////////////////////////////////////////////
///	checks whether the given line-segment (v1, v2) intersect the given box.
/**
 * vector_t has to feature a x, y and z component.
 *
 * \sa RayBoxIntersection
 */
template <class vector_t>
bool LineBoxIntersection(const vector_t& v1, const vector_t& v2,
						const vector_t& boxMin, const vector_t& boxMax);


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
bool TriangleTriangleIntersection(const MathVector<3>& p0, const MathVector<3>& p1,
								  const MathVector<3>& p2, const MathVector<3>& q0,
								  const MathVector<3>& q1, const MathVector<3>& q2,
								  MathVector<3>* ip1Out = NULL,
								  MathVector<3>* ip2Out = NULL);

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
bool TriangleBoxIntersection(const MathVector<3>& p0, const MathVector<3>& p1,
							 const MathVector<3>& p2,
							 const MathVector<3>& boxMin, const MathVector<3>& boxMax);
						
////////////////////////////////////////////////////////////////////////
///	checks whether two boxes intersect.
/**	Make sure that all parameters have the same size.*/
template <class vector_t>
bool BoxBoxIntersection(const vector_t& box1Min, const vector_t& box1Max,
						const vector_t& box2Min, const vector_t& box2Max);

////////////////////////////////////////////////////////////////////////
//	BoxBoundProbe
///	Returns true if the point lies inside or on the boundary of the box
/**	v.size() has to be the same or less than boxMin.size() and boxMax.size().
 *	If this condition is met, the method supports vectors of arbitrary dimension.*/
template <class vector_t>
bool BoxBoundProbe(const vector_t& v, const vector_t& boxMin,
					const vector_t& boxMax);

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

////////////////////////////////////////////////////////////////////////
///	Returns true if the point lies inside or on the boundary of a triangle
/**
 * This method should only be used for 2-dimensional vectors!
 */
template <class vector_t>
bool PointIsInsideTriangle(const vector_t& v, const vector_t& v0,
							const vector_t& v1, const vector_t& v2);
							
////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron
///	Returns true if the point lies inside or on the boundary of a tetrahedron
/**
 * This method does not care about the orientation of the tetrahedron.
 * This however makes it a little slower than a method that only
 * works for correctly orientated tetrahedrons.
 */
template <class vector_t>
bool PointIsInsideTetrahedron(const vector_t& v, const vector_t& v0, const vector_t& v1,
							  const vector_t& v2, const vector_t& v3);

////////////////////////////////////////////////////////////////////////
//	BinomialCoefficient
///	Returns the BinomialCoefficient
/**
 * Returns:
 * 	  n!
 * ---------
 * k! (n-k)!
 */
int BinomCoeff(int n, int k);

}//	end of namespace

////////////////////////////////////////////////
//	include implementation
#include "math_util_impl.hpp"

#endif

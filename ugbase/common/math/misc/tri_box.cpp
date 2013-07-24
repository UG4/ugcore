/* Copyright (C) Dan Ginsburg, 2000. 
 * All rights reserved worldwide.
 *
 * This software is provided "as is" without express or implied
 * warranties. You may freely copy and compile this source into
 * applications you distribute provided that the copyright text
 * below is included in the resulting source code, for example:
 * "Portions Copyright (C) Dan Ginsburg, 2000"
 */
/*
 * This Code has been modified by Sebastian Reiter.
 * Changes were made to the underlying types and to the 
 * return-value of TriangleBoxIntersection.
 * The original source-code was distributed with the book
 * 'Game Programming Gems'.
 */

////////////////////////////////////////////////////////////////////////////////////////
//
//	TriBox.cpp
//
//	Description:
//
//		This Triangle-Box intersection code was adapated from the Graphics Gems III
//		source archive 'triangleCube.c' available at:
//
//			http://www.acm.org/tog/GraphicsGems/
//
//		The main modification is that the original code performed only intersection
//		with a voxel (that is, a unit cube centered at the origin).  This code is
//		modified to take an arbitrary box and triangle and perform the scale/translation
//		necessary to get the triangle into voxel space.
//
//
#include <math.h>
#include "../ugmath.h"

///
//	Macros
//
#define INSIDE		0
#define OUTSIDE		1

#ifndef FALSE
	#define FALSE 0
#endif
#ifndef TRUE
	#define TRUE 1
#endif

#define EPS 10e-5
#define SIGN3( A ) \
	  (((A).x() < EPS) ? 4 : 0 | ((A).x() > -EPS) ? 32 : 0 | \
	   ((A).y() < EPS) ? 2 : 0 | ((A).y() > -EPS) ? 16 : 0 | \
	   ((A).z() < EPS) ? 1 : 0 | ((A).z() > -EPS) ? 8 : 0)

#define CROSS( A, B, C ) { \
  (C).x() =  (A).y() * (B).z() - (A).z() * (B).y(); \
  (C).y() = -(A).x() * (B).z() + (A).z() * (B).x(); \
  (C).z() =  (A).x() * (B).y() - (A).y() * (B).x(); \
   }
#define SUB( A, B, C ) { \
  (C).x() =  (A).x() - (B).x(); \
  (C).y() =  (A).y() - (B).y(); \
  (C).z() =  (A).z() - (B).z(); \
   }
#define LERP( A, B, C) ((B)+(A)*((C)-(B)))
#define MIN3(a,b,c) ((((a)<(b))&&((a)<(c))) ? (a) : (((b)<(c)) ? (b) : (c)))
#define MAX3(a,b,c) ((((a)>(b))&&((a)>(c))) ? (a) : (((b)>(c)) ? (b) : (c)))


namespace ug
{

//	a very simple triangle type
struct TRI
{
	vector3 m_P[3];
};

//////////////////////////////////////////////////////////////////////////////////////////
//
//	Private Functions
//
//

///
//	FacePlane()
//   
//		Which of the six face-plane(s) is point P outside of? 
//
static
int FacePlane(const vector3& p)
{
	int outcode;

	outcode = 0;
	if (p.x() >  .5) outcode |= 0x01;
	if (p.x() < -.5) outcode |= 0x02;
	if (p.y() >  .5) outcode |= 0x04;
	if (p.y() < -.5) outcode |= 0x08;
	if (p.z() >  .5) outcode |= 0x10;
	if (p.z() < -.5) outcode |= 0x20;
   
	return(outcode);
}


///
//	Bevel2d()
//
//	 Which of the twelve edge plane(s) is point P outside of? 
//
static
int Bevel2d(const vector3& p)
{
	int outcode;

	outcode = 0;
	if ( p.x() + p.y() > 1.0) outcode |= 0x001;
	if ( p.x() - p.y() > 1.0) outcode |= 0x002;
	if (-p.x() + p.y() > 1.0) outcode |= 0x004;
	if (-p.x() - p.y() > 1.0) outcode |= 0x008;
	if ( p.x() + p.z() > 1.0) outcode |= 0x010;
	if ( p.x() - p.z() > 1.0) outcode |= 0x020;
	if (-p.x() + p.z() > 1.0) outcode |= 0x040;
	if (-p.x() - p.z() > 1.0) outcode |= 0x080;
	if ( p.y() + p.z() > 1.0) outcode |= 0x100;
	if ( p.y() - p.z() > 1.0) outcode |= 0x200;
	if (-p.y() + p.z() > 1.0) outcode |= 0x400;
	if (-p.y() - p.z() > 1.0) outcode |= 0x800;
	return(outcode);
}

///
//	Bevel3d()
//	
//		Which of the eight corner plane(s) is point P outside of? 
//
static
int Bevel3d(const vector3& p)
{
	int outcode;

	outcode = 0;
	if (( p.x() + p.y() + p.z()) > 1.5) outcode |= 0x01;
	if (( p.x() + p.y() - p.z()) > 1.5) outcode |= 0x02;
	if (( p.x() - p.y() + p.z()) > 1.5) outcode |= 0x04;
	if (( p.x() - p.y() - p.z()) > 1.5) outcode |= 0x08;
	if ((-p.x() + p.y() + p.z()) > 1.5) outcode |= 0x10;
	if ((-p.x() + p.y() - p.z()) > 1.5) outcode |= 0x20;
	if ((-p.x() - p.y() + p.z()) > 1.5) outcode |= 0x40;
	if ((-p.x() - p.y() - p.z()) > 1.5) outcode |= 0x80;
	return(outcode);
}

///
//	CheckPoint()
//
//	Test the point "alpha" of the way from P1 to P2 
//  See if it is on a face of the cube   
//	Consider only faces in "mask"                   
static
int CheckPoint(const vector3& p1, const vector3& p2, number alpha, long mask)
{
	vector3 plane_point;

	plane_point.x() = LERP(alpha, p1.x(), p2.x());
	plane_point.y() = LERP(alpha, p1.y(), p2.y());
	plane_point.z() = LERP(alpha, p1.z(), p2.z());
	return(FacePlane(plane_point) & mask);
}

///
//	CheckLine()
//
//		Compute intersection of P1 --> P2 line segment with face planes 
//		Then test intersection point to see if it is on cube face       
//		Consider only face planes in "outcode_diff"                     
//		Note: Zero bits in "outcode_diff" means face line is outside of */
//
static
int CheckLine(const vector3& p1, const vector3& p2, int outcode_diff)
{

   if ((0x01 & outcode_diff) != 0)
      if (CheckPoint(p1,p2,( .5f-p1.x())/(p2.x()-p1.x()),0x3e) == INSIDE) return(INSIDE);
   if ((0x02 & outcode_diff) != 0)
      if (CheckPoint(p1,p2,(-.5f-p1.x())/(p2.x()-p1.x()),0x3d) == INSIDE) return(INSIDE);
   if ((0x04 & outcode_diff) != 0) 
      if (CheckPoint(p1,p2,( .5f-p1.y())/(p2.y()-p1.y()),0x3b) == INSIDE) return(INSIDE);
   if ((0x08 & outcode_diff) != 0) 
      if (CheckPoint(p1,p2,(-.5f-p1.y())/(p2.y()-p1.y()),0x37) == INSIDE) return(INSIDE);
   if ((0x10 & outcode_diff) != 0) 
      if (CheckPoint(p1,p2,( .5f-p1.z())/(p2.z()-p1.z()),0x2f) == INSIDE) return(INSIDE);
   if ((0x20 & outcode_diff) != 0) 
      if (CheckPoint(p1,p2,(-.5f-p1.z())/(p2.z()-p1.z()),0x1f) == INSIDE) return(INSIDE);
   return(OUTSIDE);
}

///
//	PointTriangleIntersection()
//
//		Test if 3D point is inside 3D triangle 
static
int PointTriangleIntersection(const vector3& p, const TRI& t)
{
	int		sign12,sign23,sign31;
	vector3 vect12,vect23,vect31,vect1h,vect2h,vect3h;
	vector3 cross12_1p,cross23_2p,cross31_3p;

	///
	//	First, a quick bounding-box test:                               
	//  If P is outside triangle bbox, there cannot be an intersection. 
	//
	if (p.x() > MAX3(t.m_P[0].x(), t.m_P[1].x(), t.m_P[2].x())) return(OUTSIDE);  
	if (p.y() > MAX3(t.m_P[0].y(), t.m_P[1].y(), t.m_P[2].y())) return(OUTSIDE);
	if (p.z() > MAX3(t.m_P[0].z(), t.m_P[1].z(), t.m_P[2].z())) return(OUTSIDE);
	if (p.x() < MIN3(t.m_P[0].x(), t.m_P[1].x(), t.m_P[2].x())) return(OUTSIDE);
	if (p.y() < MIN3(t.m_P[0].y(), t.m_P[1].y(), t.m_P[2].y())) return(OUTSIDE);
	if (p.z() < MIN3(t.m_P[0].z(), t.m_P[1].z(), t.m_P[2].z())) return(OUTSIDE);

	///
	//	For each triangle side, make a vector out of it by subtracting vertexes; 
	//	make another vector from one vertex to point P.                          
	//  The crossproduct of these two vectors is orthogonal to both and the      
	//  signs of its X,Y,Z components indicate whether P was to the inside or    
	//  to the outside of this triangle side.                                    
	//
	SUB(t.m_P[0], t.m_P[1], vect12);
	SUB(t.m_P[0], p,		  vect1h);	
	CROSS(vect12, vect1h, cross12_1p)
	sign12 = SIGN3(cross12_1p);      /* Extract X,Y,Z signs as 0..7 or 0...63 integer */

	SUB(t.m_P[1], t.m_P[2], vect23)
	SUB(t.m_P[1],    p, vect2h);
	CROSS(vect23, vect2h, cross23_2p)
	sign23 = SIGN3(cross23_2p);

	SUB(t.m_P[2], t.m_P[0], vect31)
	SUB(t.m_P[2],    p, vect3h);
	CROSS(vect31, vect3h, cross31_3p)
	sign31 = SIGN3(cross31_3p);

	///
	//	If all three crossproduct vectors agree in their component signs, /
	//  then the point must be inside all three.                           
	//  P cannot be OUTSIDE all three sides simultaneously.                
	//
	return (((sign12 & sign23 & sign31) == 0) ? OUTSIDE : INSIDE);
}

///
//	TriCubeIntersection()
//
//		Triangle t is compared with a unit cube,  
//		centered on the origin.                    
//		It returns INSIDE (0) or OUTSIDE(1) if t   
//		Intersects or does not intersect the cube. 
static
int TriCubeIntersection(const TRI& t)
{
	int v1_test,v2_test,v3_test;
	number d;
	vector3 vect12,vect13,norm;
	vector3 hitpp,hitpn,hitnp,hitnn;

	///
	//	First compare all three vertexes with all six face-planes 
	//	If any vertex is inside the cube, return immediately!     
	//
	if ((v1_test = FacePlane(t.m_P[0])) == INSIDE) return(INSIDE);
	if ((v2_test = FacePlane(t.m_P[1])) == INSIDE) return(INSIDE);
	if ((v3_test = FacePlane(t.m_P[2])) == INSIDE) return(INSIDE);

	///
	//	If all three vertexes were outside of one or more face-planes,
	//	return immediately with a trivial rejection!              
	//
	if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);

	///
	//	Now do the same trivial rejection test for the 12 edge planes 
	//
	v1_test |= Bevel2d(t.m_P[0]) << 8; 
	v2_test |= Bevel2d(t.m_P[1]) << 8; 
	v3_test |= Bevel2d(t.m_P[2]) << 8;
	if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);  

	///
	//	Now do the same trivial rejection test for the 8 corner planes
	//
	v1_test |= Bevel3d(t.m_P[0]) << 24; 
	v2_test |= Bevel3d(t.m_P[1]) << 24; 
	v3_test |= Bevel3d(t.m_P[2]) << 24; 
	if ((v1_test & v2_test & v3_test) != 0) return(OUTSIDE);   

	///
	//	If vertex 1 and 2, as a pair, cannot be trivially rejected 
	//  by the above tests, then see if the v1-->v2 triangle edge  
	//	intersects the cube.  Do the same for v1-->v3 and v2-->v3. 
	//  Pass to the intersection algorithm the "OR" of the outcode 
	//  bits, so that only those cube faces which are spanned by   
	//  each triangle edge need be tested.                         
	//
	if ((v1_test & v2_test) == 0)
		if (CheckLine(t.m_P[0],t.m_P[1],v1_test|v2_test) == INSIDE) return(INSIDE);
	if ((v1_test & v3_test) == 0)
		if (CheckLine(t.m_P[0],t.m_P[2],v1_test|v3_test) == INSIDE) return(INSIDE);
	if ((v2_test & v3_test) == 0)
		if (CheckLine(t.m_P[1],t.m_P[2],v2_test|v3_test) == INSIDE) return(INSIDE);

	///
	//	By now, we know that the triangle is not off to any side,     
	//	 and that its sides do not penetrate the cube.  We must now   
	//	 test for the cube intersecting the interior of the triangle. 
	//	 We do this by looking for intersections between the cube     
	//	 diagonals and the triangle...first finding the intersection  
	//	 of the four diagonals with the plane of the triangle, and    
	//	 then if that intersection is inside the cube, pursuing       
	//	 whether the intersection point is inside the triangle itself. 

	//	 To find plane of the triangle, first perform crossproduct on  
	//	 two triangle side vectors to compute the normal vector.                                       
	SUB(t.m_P[0],t.m_P[1],vect12);
	SUB(t.m_P[0],t.m_P[2],vect13);
	CROSS(vect12,vect13,norm)

	///
	//	 The normal vector "norm" X,Y,Z components are the coefficients 
	//		 of the triangles AX + BY + CZ + D = 0 plane equation.  If we   
	//	 solve the plane equation for X=Y=Z (a diagonal), we get        
	//	 -D/(A+B+C) as a metric of the distance from cube center to the 
	//	 diagonal/plane intersection.  If this is between -0.5 and 0.5, 
	//	 the intersection is inside the cube.  If so, we continue by    
	//	 doing a point/triangle intersection.                           
	//	 Do this for all four diagonals.                                
	d = norm.x() * t.m_P[0].x() + norm.y() * t.m_P[0].y() + norm.z() * t.m_P[0].z();
	hitpp.x() = hitpp.y() = hitpp.z() = d / (norm.x() + norm.y() + norm.z());
	if (fabs(hitpp.x()) <= 0.5)
		if (PointTriangleIntersection(hitpp,t) == INSIDE) return(INSIDE);
	hitpn.z() = -(hitpn.x() = hitpn.y() = d / (norm.x() + norm.y() - norm.z()));
	if (fabs(hitpn.x()) <= 0.5)
		if (PointTriangleIntersection(hitpn,t) == INSIDE) return(INSIDE);
	hitnp.y() = -(hitnp.x() = hitnp.z() = d / (norm.x() - norm.y() + norm.z()));
	if (fabs(hitnp.x()) <= 0.5)
		if (PointTriangleIntersection(hitnp,t) == INSIDE) return(INSIDE);
	hitnn.y() = hitnn.z() = -(hitnn.x() = d / (norm.x() - norm.y() - norm.z()));
	if (fabs(hitnn.x()) <= 0.5)
		if (PointTriangleIntersection(hitnn,t) == INSIDE) return(INSIDE);

	///
	//	No edge touched the cube; no cube diagonal touched the triangle. 
	//  We're done...there was no intersection.                          
	//
	return(OUTSIDE);
}

//////////////////////////////////////////////////////////////////////////////
//
//	Public Functions
//
//

///
//	TriangleBoxIntersection()
//
//		Determine if a bounding box and triangle intersect
//
bool TriangleBoxIntersection(const MathVector<3>& p0, const MathVector<3>& p1,
							 const MathVector<3>& p2,
							 const MathVector<3>& boxMin, const MathVector<3>& boxMax)
{
	vector3 Trans;
	vector3 Scale(1.0, 1.0, 1.0);
	vector3 TransMax;
	TRI		TestTri;

	///
	//	Compute the scale and transform required to make BBox
	//	a voxel
	//
	Trans.x() = (boxMax.x() + boxMin.x()) / 2;
	Trans.y() = (boxMax.y() + boxMin.y()) / 2;
	Trans.z() = (boxMax.z() + boxMin.z()) / 2;

	VecSubtract(TransMax, boxMax, Trans);
	
	if(TransMax.x() != 0)
		Scale.x() = 0.5f / TransMax.x();
	if(TransMax.y() != 0)
		Scale.y() = 0.5f / TransMax.y();
	if(TransMax.z() != 0)
		Scale.z() = 0.5f / TransMax.z();

	///
	//	Put the triangle in voxel space
	//	
	TestTri.m_P[0].x() = (p0.x() - Trans.x()) * Scale.x();
	TestTri.m_P[0].y() = (p0.y() - Trans.y()) * Scale.y();
	TestTri.m_P[0].z() = (p0.z() - Trans.z()) * Scale.z();

	TestTri.m_P[1].x() = (p1.x() - Trans.x()) * Scale.x();
	TestTri.m_P[1].y() = (p1.y() - Trans.y()) * Scale.y();
	TestTri.m_P[1].z() = (p1.z() - Trans.z()) * Scale.z();
	
	TestTri.m_P[2].x() = (p2.x() - Trans.x()) * Scale.x();
	TestTri.m_P[2].y() = (p2.y() - Trans.y()) * Scale.y();
	TestTri.m_P[2].z() = (p2.z() - Trans.z()) * Scale.z();

	///
	//	Test against the voxel
	//
	return(TriCubeIntersection(TestTri) == INSIDE);	
}

}//	end of namespace

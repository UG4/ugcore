//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d14

#ifndef __H__UGMATH__MATH_UTIL_IMPL__
#define __H__UGMATH__MATH_UTIL_IMPL__

#include <algorithm>
#include "common/common.h"
#include "common/math/math_vector_matrix/math_vector_functions.h"
#include "common/math/math_vector_matrix/math_matrix_functions.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	deg_to_rad
template <class TNumber>
inline TNumber
deg_to_rad(TNumber deg)
{
	return deg * PI / 180.;
}

////////////////////////////////////////////////////////////////////////
//	rad_to_deg
template <class TNumber>
inline TNumber
rad_to_deg(TNumber rad)
{
	return rad * 180. / PI;
}

////////////////////////////////////////////////////////////////////////
//	urand
template <class TNumber>
TNumber
urand(TNumber lowerBound, TNumber upperBound)
{
	long t = rand();
	if(t == RAND_MAX)
		t -= 1;

	return lowerBound + (TNumber)((upperBound - lowerBound) * ((double)t / (double)RAND_MAX));
}

////////////////////////////////////////////////////////////////////////
//	clip
template <class TNumber>
TNumber
clip(TNumber val, TNumber lowerBound, TNumber upperBound)
{
	if(val > upperBound)
		return upperBound;
	else if(val < lowerBound)
		return lowerBound;
	return val;
}

////////////////////////////////////////////////////////////////////////
template <class TNumber>
inline TNumber sq(TNumber val)
{
	return val * val;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
void CalculateCenter(vector_t& centerOut, const vector_t* pointSet,
					 size_t numPoints)
{
	for(size_t i = 0; i < centerOut.size(); ++i)
		centerOut[i] = 0;
	
	if(numPoints > 0){
		for(size_t i = 0; i < numPoints; ++i)
			VecAdd(centerOut, centerOut, pointSet[i]);
	
		VecScale(centerOut, centerOut, 1. / (number)numPoints);
	}
}

template <class vector_t>
vector_t
TriangleBarycenter(const vector_t& p1, const vector_t& p2, const vector_t& p3)
{
	vector_t bc;
	VecScaleAdd(bc, 1./3., p1, 1./3., p2, 1./3., p3);
	return bc;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number DropAPerpendicular(vector_t& vOut, const vector_t& v,
							const vector_t& v0, const vector_t& v1)
{
//	project v onto v' on the edge (v0, v1) so that (v'-v)*(v0-v1) = 0
	vector_t e[2];
	VecSubtract(e[0], v, v0);
	VecSubtract(e[1], v1, v0);

	number d1 = VecDot(e[0], e[1]);
	number d2 = VecDot(e[1], e[1]);
	
//	avoid division by zero
	if(fabs(d2) > SMALL)
	{
	//	calculate the projection p'
		number s = d1/d2;
		VecScale(e[1], e[1], s);
		VecAdd(vOut, v0, e[1]);
		return s;
	}
	else
		vOut = v0;
	return 0;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number ProjectPointToRay(vector_t& vOut, const vector_t& v,
							const vector_t& from, const vector_t& dir)
{
	vector_t tmpDir;
	VecSubtract(tmpDir, v, from);

	number d1 = VecDot(tmpDir, dir);
	number d2 = VecDot(dir, dir);
	
//	avoid division by zero
	if(fabs(d2) > SMALL)
	{
	//	calculate the projection p'
		number s = d1/d2;
		VecScale(tmpDir, dir, s);
		VecAdd(vOut, from, tmpDir);
		return s;
	}
	else{
		vOut = from;
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number ProjectPointToLine(vector_t& vOut, const vector_t& v,
						  const vector_t& from, const vector_t& to)
{
	vector_t dir;
	VecSubtract(dir, to, from);
	return ProjectPointToRay(vOut, v, from, dir);
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
inline
number DistancePointToLine(const vector_t& v, const vector_t& v1,
						  const vector_t& v2)
{
	number t;
	return DistancePointToLine(t, v, v1, v2);
}

template <class vector_t>
number DistancePointToLine(number& tOut, const vector_t& v,
						   const vector_t& v1, const vector_t& v2)
{
	vector_t tmp;
	tOut = DropAPerpendicular(tmp, v, v1, v2);
	if(tOut > 1){
		tOut = 1;
		return VecDistance(v, v2);
	}
	else if(tOut < 0){
		tOut = 0;
		return VecDistance(v, v1);
	}
	else
		return VecDistance(v, tmp);
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
inline
number DistancePointToRay(const vector_t& v, const vector_t& from,
						  const vector_t& dir)
{
	vector_t tmp;
	ProjectPointToRay(tmp, v, from, dir);
	return VecDistance(v, tmp);
}

template <class vector_t>
inline
number DistancePointToRay(vector_t& vOut, number& tOut, const vector_t& v,
						  const vector_t& from, const vector_t& dir)
{
	tOut = ProjectPointToRay(vOut, v, from, dir);
	return VecDistance(v, vOut);
}

template <class vector_t>
number DistancePointToTriangle(vector_t& vOut, number& bc1Out, number& bc2Out,
							const vector_t& p, const vector_t& v1, const vector_t& v2,
							const vector_t& v3, const vector_t& n)
{
//	parameter of line-intersection and of point-ray-projection.
	number t;
	
//	first try an orthogonal projection of p onto the triangle
	if(RayTriangleIntersection(vOut, bc1Out, bc2Out, t, v1, v2, v3, p, n))
	{	//min distance found
		return VecDistance(vOut, p);
	}
	
//	if ortho projection is not in triangle perform edge-point distance
	number d, tmpDist, tmpT;
	vector_t vDir, vTmp;
	int bestIndex = 0;
	
	VecSubtract(vDir, v2, v1);
	d = DistancePointToRay(vOut, t, p, v1, vDir);
	bc1Out = t; bc2Out = 0;
	
	VecSubtract(vDir, v3, v1);
	tmpDist = DistancePointToRay(vTmp, tmpT, p, v1, vDir);
	if(tmpDist < d){
		bestIndex = 1;
		d = tmpDist;
		t = tmpT;
		bc1Out = 0; bc2Out = tmpT;
		vOut = vTmp;
	}
	
	VecSubtract(vDir, v3, v2);
	tmpDist = DistancePointToRay(vTmp, tmpT, p, v2, vDir);
	if(tmpDist < d){
		bestIndex = 2;
		d = tmpDist;
		t = tmpT;
		bc1Out = 1. - t; bc2Out = t;
		vOut = vTmp;
	}
	
//	we now have to check whether the projection to an edge was
//	orthogonal. If not we'll return the distance to a point.
	if(t > 0 && t < 1){
		return d;
	}
	else{
		switch(bestIndex){
			case 0:
				if(t < 0.5){
					vOut = v1;
					bc1Out = 0; bc2Out = 0;
				}
				else{
					vOut = v2;
					bc1Out = 1; bc2Out = 0;
				}
				break;
				
			case 1:
				if(t < 0.5){
					vOut = v1;
					bc1Out = 0; bc2Out = 0;
				}
				else{
					vOut = v3;
					bc1Out = 0; bc2Out = 1;
				}
				break;
			case 2:
				if(t < 0.5){
					vOut = v2;
					bc1Out = 1; bc2Out = 0;
				}
				else{
					vOut = v3;
					bc1Out = 0; bc2Out = 1;
				}
				break;
		}
	}
	
//	return the distance
	return VecDistance(p, vOut);
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number DistancePointToPlane(const vector_t& v, const vector_t& p,
							const vector_t& n)
{
	vector_t vTmp;
	ProjectPointToPlane(vTmp, v, p, n);
	return VecDistance(vTmp, v);
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
void ProjectPointToPlane(vector_t& vOut, const vector_t& v,
						const vector_t& p, const vector_t& n)
{
//	the vector from p to v
	vector_t t, norm;
	VecSubtract(t, v, p);

	VecNormalize(norm, n);

//	scale the normal with the dot-product with the direction
	VecScale(t, norm, VecDot(norm, t));

//	subtract the scaled normal from the original point
	VecSubtract(vOut, v, t);
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool RayPlaneIntersection(vector_t& vOut, number& tOut,
						  const vector_t& rayFrom, const vector_t& rayDir,
						  const vector_t& p, const vector_t& n)
{
//	solve: t = (p-rayFrom)*n / rayDir*n
	number denom = VecDot(rayDir, n);
	if(fabs(denom) < SMALL_SQ)
		return false;
	
//	calculate intersection parameter
	vector_t v;
	VecSubtract(v, p, rayFrom);
	tOut = VecDot(v, n) / denom;
	
//	calculate intersection point
	VecScale(v, rayDir, tOut);
	VecAdd(vOut, rayFrom, v);
	return true;
}


////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool RayRayIntersection2d(vector_t &vOut, number& t0Out, number& t1Out,
						   const vector_t &p0, const vector_t &dir0,
						   const vector_t &p1, const vector_t &dir1)
{
	// we search for the intersection of the ray vFrom + tOut*vDir with the
	// Line p0 + bcOut*(p1-p0). Intersection is true, if c0 in [0,1].
	// We set up the system
	//   | dir0[0] , - dir1[0] | ( t0Out ) = ( (p1 - p0)[0] )
	//   | dir0[1] , - dir1[1] | ( t1Out ) = ( (p1 - p0)[1] )

	vector_t v, b;
	MathMatrix<2,2> m, mInv;

	// set up matrix
	m[0][0] = dir0[0]; 	m[0][1] = -1 * dir1[0];
	m[1][0] = dir0[1]; 	m[1][1] = -1 * dir1[1];

	// invert matrix
	number det = Determinant(m);

	// if det == 0.0 lines are parallel
	if(det == 0.0) return false;

	// compute inverse
	Inverse(mInv, m);

	// compute rhs of system
	VecSubtract(b, p1, p0);

	// solve system
	MatVecMult(v, mInv, b);

	// parameters of the intersection
	t0Out = v[0];
	t1Out = v[1];

	// compute intersection point
	vOut = p0;
	VecScaleAppend(vOut, t0Out, dir0);

	return true;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool RayLineIntersection2d(vector_t &vOut, number& bcOut, number& tOut,
						   const vector_t &p0, const vector_t &p1,
						   const vector_t &vFrom, const vector_t &vDir)
{
	vector_t dir0;
	VecSubtract(dir0, p1, p0);
	if(RayRayIntersection2d(vOut, bcOut, tOut, p0, dir0, vFrom, vDir)){
		if(bcOut >= 0 && bcOut <= 1.)
			return true;
	}
	return false;
}

template <class vector_t>
bool LineLineIntersection2d(vector_t &vOut, number& t0Out, number& t1Out,
						   const vector_t &from0, const vector_t &to0,
						   const vector_t &from1, const vector_t &to1,
						   const number threshold)
{
	vector_t dir0, dir1;
	VecSubtract(dir0, to0, from0);
	VecSubtract(dir1, to1, from1);
	if(RayRayIntersection2d(vOut, t0Out, t1Out, from0, dir0, from1, dir1)){
		if((t0Out >= -threshold) && (t0Out <= (1. + threshold))
			&& (t1Out >= -threshold) && (t1Out <= (1. + threshold)))
		{
			return true;
		}
	}
	return false;
}

template <class vector_t>
bool RayRayProjection(number& t1Out, number& t2Out,
						const vector_t& from1, const vector_t& dir1,
						const vector_t& from2, const vector_t& dir2)
{
	vector_t ab;
	VecSubtract(ab, from2, from1);
	number l11 = VecDot(dir1, dir1);
	number l12 = -VecDot(dir1, dir2);
	number l22 = VecDot(dir2, dir2);
	number ra = VecDot(dir1, ab);
	number rb = -VecDot(dir2, ab);

//	l11 and l22 are always >= 0
	if((l11 < SMALL_SQ) || (l22 < SMALL_SQ))
		return false;

	number tmp = l11 * l22 - l12 * l12;
	if(fabs(tmp) < SMALL)
		return false;

	t2Out = (l11*rb - l12*ra) / tmp;
	t1Out = (ra - l12*t2Out) / l11;
	return true;
}

template <class vector_t>
bool LineLineProjection(number& t1Out, number& t2Out,
						  const vector_t& a1, const vector_t& a2,
						  const vector_t& b1, const vector_t& b2)
{
	vector_t dirA, dirB;
	VecSubtract(dirA, a2, a1);
	VecSubtract(dirB, b2, b1);

	if(RayRayProjection(t1Out, t2Out, a1, dirA, b1, dirB)){
		if((t1Out >= 0) && (t1Out <= 1.) && (t2Out >= 0) && (t2Out <= 1.))
			return true;
		return false;
	}
	return false;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool RayTriangleIntersection(vector_t &vOut, number& bc1Out, number& bc2Out, number& tOut,
						   const vector_t &p0, const vector_t &p1, const vector_t &p2, 
						   const vector_t &vFrom, const vector_t &vDir)
{
//	this piece of code is rather old and should probably be replaced by a
//	new method to improve speed and robustness.
//	The current implementation solves a 3x3 linear system using gaussian
//	elimination with pivoting to find the intersection of the ray with
//	the plane in which the triangle lies.
//	The local coordinates r and s, that describe the position of the
//	intersection point in the planes parameter-form, are then examined to
//	check whether the intersection-point lies inside the triangle or not
//	(r and s can at the same time be seen as the barycentric coordinates of
//	the triangle).
//	Division by zero is avoided (currently a == check is used - should be 
//	replaced by a comparision to SMALL)

/*
	p2
	|\
  v1|  \v2
	|____\
  p0  v0  p1


	r*v0 + s*v1 - t* vDir = b

	m00	m01	m02
	m10	m11	m12
	m20	m21	m22
*/

	number m[3][3];
	number b[3];

	m[0][0] = p1.x() - p0.x();	m[0][1] = p2.x() - p0.x();	m[0][2] = -vDir.x();	b[0] = vFrom.x() - p0.x();
	m[1][0] = p1.y() - p0.y();	m[1][1] = p2.y() - p0.y();	m[1][2] = -vDir.y();	b[1] = vFrom.y() - p0.y();
	m[2][0] = p1.z() - p0.z();	m[2][1] = p2.z() - p0.z();	m[2][2] = -vDir.z();	b[2] = vFrom.z() - p0.z();

	int i1, i2, i3, j;
	number fac;
	bc1Out = 0;
	bc2Out = 0;
	tOut = 0;

	number dBestEntry = fabs(m[0][0]);
	i1 = 0;
	fac = fabs(m[1][0]);
	if(fac > dBestEntry){
		dBestEntry = fac;
		i1 = 1;
	}
	
	if(fabs(m[2][0]) > dBestEntry)
		i1 = 2;


	if(m[i1][0]){
		for(i2 = 0; i2 < 3; ++i2){
			if(i2!=i1){
				if(m[i2][0]){
					fac = -m[i2][0]/m[i1][0];
					for(j = 0; j < 3; ++j)
						m[i2][j] = m[i2][j] + fac*m[i1][j];
					b[i2] = b[i2] + fac * b[i1];
				}
			}
		}
		
		i2 = (i1 + 1) % 3;
		i3 = (i1 + 2) % 3;
		if(fabs(m[i2][1]) < fabs(m[i3][1])){
			int ti = i2;
			i2 = i3;
			i3 = ti;
		}
		
		if((m[i2][1]!=0) && (m[i3][1]!=0)){
			fac = -m[i3][1]/m[i2][1];
			for(j = 1; j < 3; ++j)
				m[i3][j] = m[i3][j] + fac*m[i2][j];
			b[i3] = b[i3] + fac * b[i2];
		}
		
		//calculate tOut (t)
		if(m[i3][2])
			tOut = b[i3] / m[i3][2];
		else if(b[i3] != 0)
			return false;

		//calculate bc2Out (s)
		b[i2] -= tOut*m[i2][2];
		if(m[i2][1])
			bc2Out = b[i2] / m[i2][1];
		else if(b[i2] != 0)
			return false;

		//calculate bc1Out (r)
		b[i1] -= (tOut*m[i1][2] + bc2Out*m[i1][1]);
		if(m[i1][0])
			bc1Out = b[i1] / m[i1][0];
		else if(b[i1] != 0)
			return false;

		if((bc1Out >=-SMALL ) && (bc2Out >= -SMALL ) && ((bc1Out + bc2Out) <= 1.+SMALL))
		{
			vOut.x() = vFrom.x() + tOut*vDir.x();
			vOut.y() = vFrom.y() + tOut*vDir.y();
			vOut.z() = vFrom.z() + tOut*vDir.z();
			return true;
		}
	}
	return false;
}

template <class vector_t> inline
bool RayTriangleIntersection(vector_t &vOut, const vector_t &p0,
						   const vector_t &p1, const vector_t &p2, 
						   const vector_t &vFrom, const vector_t &vDir)
{
	number r, s, t;
	return RayTriangleIntersection(vOut, r, s, t, p0, p1, p2, vFrom, vDir);
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool RayBoxIntersection(const vector_t& rayFrom, const vector_t& rayDir,
						const vector_t& boxMin, const vector_t& boxMax,
						number* tNearOut, number* tFarOut)
{
	number tMin = 0, tMax = 0;// initialized only to avoid mislead compiler warnings...
	number t1, t2;
	bool bMinMaxSet = false;
	
	if(fabs(rayDir.x()) > SMALL)
	{
		//get xNear and xFar
		t1 = (boxMin.x() - rayFrom.x()) / rayDir.x();
		t2 = (boxMax.x() - rayFrom.x()) / rayDir.x();
		if(t1 > t2)
			std::swap(t1, t2);
		tMin	= t1;
		tMax	= t2;
		bMinMaxSet = true;
	}
	else
	{
		if(rayFrom.x() < boxMin.x())
			return false;
		if(rayFrom.x() > boxMax.x())
			return false;
	}
	
	if(fabs(rayDir.y()) > SMALL)
	{
		//get yNear and yFar
		t1 = (boxMin.y() - rayFrom.y()) / rayDir.y();
		t2 = (boxMax.y() - rayFrom.y()) / rayDir.y();
		if(t1 > t2)
			std::swap(t1, t2);
		if(bMinMaxSet)
		{
			if((t1 <= tMax) && (t2 >= tMin))
			{
				tMin = std::max(t1, tMin);
				tMax = std::min(t2, tMax);
			}
			else
				return false;
		}
		else
		{
			tMin	= t1;
			tMax	= t2;
		}
		bMinMaxSet = true;
	}
	else
	{
		if(rayFrom.y() < boxMin.y())
			return false;
		if(rayFrom.y() > boxMax.y())
			return false;		
	}
	
	if(fabs(rayDir.z()) > SMALL)
	{
		//get zNear and zFar
		t1 = (boxMin.z() - rayFrom.z()) / rayDir.z();
		t2 = (boxMax.z() - rayFrom.z()) / rayDir.z();
		if(t1 > t2)
			std::swap(t1, t2);
		if(bMinMaxSet)
		{
			if((t1 <= tMax) && (t2 >= tMin))
			{
				tMin = std::max(t1, tMin);
				tMax = std::min(t2, tMax);
			}
			else
				return false;
		}
		else
		{
			tMin	= t1;
			tMax	= t2;
		}
		bMinMaxSet = true;
	}
	else
	{
		if(rayFrom.z() < boxMin.z())
			return false;
		if(rayFrom.z() > boxMax.z())
			return false;		
	}
	
	if(bMinMaxSet)
	{
	//	the ray intersects the box
	//	lets calculate tNear and tFar and return true
		if(fabs(tMin) > fabs(tMax))
			std::swap(tMin, tMax);
		if(tNearOut)
			*tNearOut = tMin;
		if(tFarOut)
			*tFarOut = tMax;
		return true;
	}
	else
	{
	//	the ray has no direction -> we'll check if the From-point lies inside the box
		if(BoxBoundProbe(rayFrom, boxMin, boxMax)){
		
			if(*tNearOut)
				*tNearOut = 0;
			if(*tFarOut)
				*tFarOut = 0;
			return true;
		}
		else
			return false;
	}
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool LineBoxIntersection(const vector_t& v1, const vector_t& v2,
						const vector_t& boxMin, const vector_t& boxMax)
{
	number tNear, tFar;

	vector_t vDir;
	VecSubtract(vDir, v2, v1);
	if(RayBoxIntersection(v1, vDir, boxMin, boxMax, &tNear, &tFar))
	{		
		return ((tNear * tFar < 0) || (tNear >= 0 && tNear <= 1.));
	}

	return false;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
int RaySphereIntersection(number& s1Out, number& s2Out,
						  const vector_t& v, const vector_t& dir,
						  const vector_t& center, number radius)
{
	vector_t p;
	number s = ProjectPointToRay(p, center, v, dir);
	number h = VecDistance(p, center);
	if(h > radius)
		return 0;
	if(h > radius - SMALL){
		s1Out = s;
		s2Out = s;
		return 1;
	}

	number dirLen = VecLength(dir);
	if(dirLen == 0){
		s1Out = s2Out = 0;
		return 0;
	}

	number a = sqrt(radius * radius - h * h);

	number sa = a / dirLen;
	s1Out = s + sa;
	s2Out = s - sa;
	return 2;
}

template <class vector_t>
int LineSphereIntersection(number& s1Out, number& s2Out,
						  const vector_t& v1, const vector_t& v2,
						  const vector_t& center, number radius)
{
	vector_t dir;
	VecSubtract(dir, v2, v1);
	int num = RaySphereIntersection(s1Out, s2Out, v1, dir, center, radius);

	switch(num){
		case 0: return 0;
		case 1: if((s1Out >= 0) && (s1Out <= 1))
					return 1;
				return 0;
		case 2:{
			if((s1Out < 0) || (s1Out > 1))
				--num;
			if((s2Out < 0) || (s2Out > 1))
				--num;
			else if(num == 1)
				s1Out = s2Out;
			return num;
		}
	}
	return 0;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
bool BoxBoxIntersection(const vector_t& box1Min, const vector_t& box1Max,
						const vector_t& box2Min, const vector_t& box2Max)
{
	for(size_t i = 0; i < box1Min.size(); ++i){
		if(box1Min[i] > box2Max[i] || box1Max[i] < box2Min[i])
			return false;
	}

	return true;
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number TriangleArea(const vector_t& p1, const vector_t& p2, const vector_t& p3)
{
//	the projection of p3 onto the line defined by p1 and p2
	vector_t v;
	DropAPerpendicular(v, p3, p1, p2);
//	calculate the area
	return 0.5 * sqrt(VecDistanceSq(p1, p2) * VecDistanceSq(v, p3));
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number GeometricApproximationDegree(vector_t& n1, vector_t& n2, vector_t& n3,
									vector_t& tn)
{
	return std::min(VecDot(n1, tn), std::min(VecDot(n2, tn), VecDot(n3, tn)));
}

////////////////////////////////////////////////////////////////////////
template <class vector_t>
number TriangleQuality_Area(const vector_t& p1, const vector_t& p2,
							const vector_t& p3)
{
	number edgeSum = VecDistanceSq(p1, p2) + VecDistanceSq(p2, p3)
					+ VecDistanceSq(p1, p3);
	if(edgeSum > SMALL)
	{
	//	4*sqrt(3) = 6.9282032
		return 6.9282032 * TriangleArea(p1, p2, p3) / edgeSum;
	}

//	a triangle whose sides have all zero length is considered to be a bad triangle.
	return 0;
}

////////////////////////////////////////////////////////////////////////
//	BoxBoundProbe
template <class vector_t>
bool BoxBoundProbe(const vector_t& v, const vector_t& boxMin,
					const vector_t& boxMax)
{
	for(size_t i = 0; i < v.size(); ++i){
		if(v[i] < boxMin[i] || v[i] > boxMax[i])
			return false;
	}
	return true;
}

////////////////////////////////////////////////////////////////////////
//	PointIsInsideTriangle
template <class vector_t>
bool PointIsInsideTriangle(const vector_t& v, const vector_t& v0,
							const vector_t& v1, const vector_t& v2)
{
//	we'll check for each side of the tri, whether v and the point of
//	the tri, that does not lie on the edge, do lie on the same side.
	vector_t e;			// the examined edge
	vector_t edgeNorm;	// the normal of the examined edge
	vector_t tv1, tv2;	// the direction of a tri-point to v
	

	VecSubtract(e, v1, v0);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v2, v0);
	VecSubtract(tv2, v, v0);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;

	VecSubtract(e, v2, v1);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v0, v1);
	VecSubtract(tv2, v, v1);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;

	VecSubtract(e, v0, v2);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v1, v2);
	VecSubtract(tv2, v, v2);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;
	
//	all tests succeeded. return true.
	return true;
}


////////////////////////////////////////////////////////////////////////
//	PointIsInsideTriangle_HighAcc
template <class vector_t>
bool PointIsInsideTriangle_HighAcc(const vector_t& v, const vector_t& v0,
								   const vector_t& v1, const vector_t& v2)
{
	using std::max;

//	we'll check for each side of the tri, whether v and the point of
//	the tri, that does not lie on the edge, do lie on the same side.
	vector_t e;			// the examined edge
	vector_t edgeNorm;	// the normal of the examined edge
	vector_t tv1, tv2;	// the direction of a tri-point to v
	typedef typename vector_t::value_type value_t;

	const value_t locSmallSq = SMALL_SQ * sqrt(max(VecDistanceSq(v0, v1),
								  	 	 	 	   max(VecDistanceSq(v1, v2),
								  	 	 	 		   VecDistanceSq(v2, v0))));

	VecSubtract(e, v1, v0);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v2, v0);
	VecSubtract(tv2, v, v0);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -locSmallSq)
		return false;

	VecSubtract(e, v2, v1);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v0, v1);
	VecSubtract(tv2, v, v1);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -locSmallSq)
		return false;

	VecSubtract(e, v0, v2);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v1, v2);
	VecSubtract(tv2, v, v2);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -locSmallSq)
		return false;

//	all tests succeeded. return true.
	return true;
}

////////////////////////////////////////////////////////////////////////
//	PointIsInsideQuadrilateral
template <class vector_t>
bool PointIsInsideQuadrilateral(const vector_t& v, const vector_t& v0,
								const vector_t& v1, const vector_t& v2,
								const vector_t& v3)
{
//	we'll check for each side of the quad, whether v and a point of
//	the quad, that does not lie on the edge, do lie on the same side.
	vector_t e;			// the examined edge
	vector_t edgeNorm;	// the normal of the examined edge
	vector_t tv1, tv2;	// the direction of a quad-point to v


	VecSubtract(e, v1, v0);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v2, v0);
	VecSubtract(tv2, v, v0);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;

	VecSubtract(e, v2, v1);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v3, v1);
	VecSubtract(tv2, v, v1);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;

	VecSubtract(e, v3, v2);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v0, v2);
	VecSubtract(tv2, v, v2);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;

	VecSubtract(e, v0, v3);
	edgeNorm.x() = e.y();
	edgeNorm.y() = -e.x();
	VecSubtract(tv1, v1, v3);
	VecSubtract(tv2, v, v3);
	if(VecDot(tv1, edgeNorm) * VecDot(tv2, edgeNorm) < -SMALL_SQ)
		return false;

//	all tests succeeded. return true.
	return true;
}

////////////////////////////////////////////////////////////////////////
//	PointIsInsideTetrahedron
template <class vector_t>
bool PointIsInsideTetrahedron(const vector_t& v, const vector_t& v0, const vector_t& v1,
							  const vector_t& v2, const vector_t& v3)
{
//	we'll check for each side of the tet, whether v and the point of
//	the tet, that does not lie in the plane, lie on the same side.
	vector_t n;			// the normal of the examined face
	vector_t e1, e2;	// directions of two edges of a face
	number pn;			// dot product of a point in the plane with the normal

//	check side 0, 2, 1
	VecSubtract(e1, v2, v0);
	VecSubtract(e2, v1, v0);
	VecCross(n, e1, e2);
	pn = VecDot(v0, n);
	if((VecDot(v3, n) - pn) * (VecDot(v, n) - pn) < -SMALL_SQ)
		return false;

//	check side 0, 1, 3
	VecSubtract(e1, v1, v0);
	VecSubtract(e2, v3, v0);
	VecCross(n, e1, e2);
	pn = VecDot(v0, n);
	if((VecDot(v2, n) - pn) * (VecDot(v, n) - pn) < -SMALL_SQ)
		return false;

//	check side 1, 2, 3
	VecSubtract(e1, v2, v1);
	VecSubtract(e2, v3, v1);
	VecCross(n, e1, e2);
	pn = VecDot(v1, n);
	if((VecDot(v0, n) - pn) * (VecDot(v, n) - pn) < -SMALL_SQ)
		return false;

//	check side 0, 3, 2
	VecSubtract(e1, v3, v0);
	VecSubtract(e2, v2, v0);
	VecCross(n, e1, e2);
	pn = VecDot(v0, n);
	if((VecDot(v1, n) - pn) * (VecDot(v, n) - pn) < -SMALL_SQ)
		return false;

//	all tests succeeded. return true.
	return true;
}

////////////////////////////////////////////////////////////////////////
//	ReflectVectorAtPlane
template <class vector_t>
void ReflectVectorAtPlane(vector_t& vReflectedOut, const vector_t& v,
                          const vector_t& n, const vector_t& r0)
{
	const number s = 2 * (VecDot(v, n) - VecDot(n, r0)) / VecDot(n, n);
	VecScaleAdd(vReflectedOut, 1.0, v, -s, n);
}

}//	end of namespace

#endif

//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m10 d14

#ifndef __H__UGMATH__MATH_UTIL_IMPL__
#define __H__UGMATH__MATH_UTIL_IMPL__

#include <algorithm>

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	urand
template <class TNumber>
TNumber
urand(TNumber lowerBound, TNumber upperBound)
{
	long t = rand();
	if(t == RAND_MAX)
		t -= 1;

	return lowerBound + (TNumber)((upperBound - lowerBound) * ((float)t / (float)RAND_MAX));
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
template <class vector_t>
number DropAPerpendicular(vector_t& vOut, const vector_t& v0,
							const vector_t& v1, const vector_t& v)
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
void ProjectPointToPlane(vector_t& vOut, const vector_t& v,
						const vector_t& p, const vector_t& n)
{
//	the vector from p to v
	vector_t t;
	VecSubtract(t, v, p);

//	scale the normal with the dot-product with the direction
	VecScale(t, n, VecDot(n, t));

//	subtract the scaled normal from the original point
	VecSubtract(vOut, v, t);
}

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

	m[0][0] = p1.x - p0.x;	m[0][1] = p2.x - p0.x;	m[0][2] = -vDir.x;	b[0] = vFrom.x - p0.x;
	m[1][0] = p1.y - p0.y;	m[1][1] = p2.y - p0.y;	m[1][2] = -vDir.y;	b[1] = vFrom.y - p0.y;
	m[2][0] = p1.z - p0.z;	m[2][1] = p2.z - p0.z;	m[2][2] = -vDir.z;	b[2] = vFrom.z - p0.z;

	int i1, i2, i3, j;
	number fac;
	number r, s, t;

	number dBestEntry = fabs(m[0][0]);
	i1 = 0;
	fac = fabs(m[1][0]);
	if(fac > dBestEntry)
	{
		dBestEntry = fac;
		i1 = 1;
	}
	if(fabs(m[2][0]) > dBestEntry)
		i1 = 2;


	if(m[i1][0])
	{
		for(i2 = 0; i2 < 3; ++i2)
		{
			if(i2!=i1)
			{
				if(m[i2][0])
				{
					fac = -m[i2][0]/m[i1][0];
					for(j = 0; j < 3; ++j)
						m[i2][j] = m[i2][j] + fac*m[i1][j];
					b[i2] = b[i2] + fac * b[i1];
				}
			}
		}
		
		i2 = (i1 + 1) % 3;
		i3 = (i1 + 2) % 3;
		if(fabs(m[i2][1]) < fabs(m[i3][1]))
		{
			int ti = i2;
			i2 = i3;
			i3 = ti;
		}
		if((m[i2][1]!=0) && (m[i3][1]!=0))
		{
			fac = -m[i3][1]/m[i2][1];
			for(j = 1; j < 3; ++j)
				m[i3][j] = m[i3][j] + fac*m[i2][j];
			b[i3] = b[i3] + fac * b[i2];
		}
		//calculate t
		if(m[i3][2])
			t = b[i3] / m[i3][2];
		else
		{
			if(b[i3] == 0)
				t = 0;
			else
				return false;
		}
		//calculate s
		b[i2] -= t*m[i2][2];
		if(m[i2][1])
			s = b[i2] / m[i2][1];
		else
		{
			if(b[i2] == 0)
				s = 0;
			else
				return false;
		}
		//calculate r
		b[i1] -= (t*m[i1][2] + s*m[i1][1]);
		if(m[i1][0])
			r = b[i1] / m[i1][0];
		else
		{
			if(b[i1] == 0)
				r = 0;
			else
				return false;
		}
	}

	if((r >=-SMALL ) && (s >= -SMALL ) && ((r + s) <= 1+SMALL))
	{
		vOut.x = vFrom.x + t*vDir.x;
		vOut.y = vFrom.y + t*vDir.y;
		vOut.z = vFrom.z + t*vDir.z;
		bc1Out = r;
		bc2Out = s;
		tOut = t;
		return true;
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
number TriangleArea(const vector_t& p1, const vector_t& p2, const vector_t& p3)
{
//	the projection of p3 onto the line defined by p1 and p2
	vector_t v;
	DropAPerpendicular(v, p1, p2, p3);
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

}//	end of namespace

#endif

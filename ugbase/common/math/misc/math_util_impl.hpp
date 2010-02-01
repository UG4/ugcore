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

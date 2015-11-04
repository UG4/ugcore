#include <vector>
#include "math_util.h"
#include "../ugmath.h"
#include "../ugmath_types.h"
#include "lineintersect_utils.h"


using namespace std;

namespace ug
{

bool TriangleCircumcenter(vector2& centerOut, const vector2& p1,
						  const vector2& p2, const vector2& p3)
{
	using std::swap;

	number d12 = VecDistanceSq(p1, p2);
	number d23 = VecDistanceSq(p2, p3);
	number d13 = VecDistanceSq(p1, p3);
/*
//	if any of the sides is too short, then we have to abort
	if(d12 < SMALL || d23 < SMALL || d13 < SMALL){
		return false;
	}
*/
//	centers of sides and side-normals
	vector2 c1, c2, n1, n2;

//	for maximal accuracy, we'll choose the two longest sides
	if(d12 >= d23){
		VecScaleAdd(c1, 0.5, p1, 0.5, p2);
		VecSubtract(n1, p2, p1);//	perform swapping later
		if(d23 >= d13){
			VecScaleAdd(c2, 0.5, p2, 0.5, p3);
			VecSubtract(n2, p3, p2);//	perform swapping later
		}
		else{
			VecScaleAdd(c2, 0.5, p1, 0.5, p3);
			VecSubtract(n2, p3, p1);//	perform swapping later
		}
	}
	else{
		VecScaleAdd(c1, 0.5, p2, 0.5, p3);
		VecSubtract(n1, p3, p2);//	perform swapping later
		if(d12 >= d13){
			VecScaleAdd(c2, 0.5, p1, 0.5, p2);
			VecSubtract(n2, p2, p1);//	perform swapping later
		}
		else{
			VecScaleAdd(c2, 0.5, p1, 0.5, p3);
			VecSubtract(n2, p3, p1);//	perform swapping later
		}
	}

//	swap normal-coefficients
	swap(n1.x(), n1.y());
	n1.x() *= -1.;
	swap(n2.x(), n2.y());
	n2.x() *= -1.;

//	calculate intersection of the two lines
	number t0, t1;
	if(!RayRayIntersection2d(centerOut, t0, t1, c1, n1, c2, n2)){
	//	line-line intersection failed. return the barycenter instead.
		return false;
	}

	return true;
}

bool TriangleCircumcenter(vector3& centerOut, const vector3& p1,
						  const vector3& p2, const vector3& p3)
{
	number d12 = VecDistanceSq(p1, p2);
	number d13 = VecDistanceSq(p1, p3);
	number d23 = VecDistanceSq(p2, p3);

//	sort p1, p2 and p3 into v1, v2, v3, so that v1 is the vertex at which the
//	two shorter edges meet
	vector3 v1, v2, v3;
	if(d12 < d13){
		if(d13 < d23){
		//	d23 is biggest
			v1 = p1; v2 = p2; v3 = p3;
		}
		else{
		//	d13 is biggest
			v1 = p2; v2 = p3; v3 = p1;
		}
	}
	else{
		if(d12 < d23){
		//	d23 is biggest
			v1 = p1; v2 = p2; v3 = p3;
		}
		else{
		//	d12 is biggest
			v1 = p3; v2 = p1; v3 = p2;
		}
	}

	vector3 dir12, dir13, dir23;
	VecSubtract(dir12, v2, v1);
	VecSubtract(dir13, v3, v1);
	VecSubtract(dir23, v3, v2);

//	we'll construct a line through the center of v1, v2, with a direction
//	normal to dir12, which lies in the triangles plane.
	vector3 proj;
	ProjectPointToRay(proj, v1, v2, dir23);
	VecSubtract(proj, proj, v1);

	number a = VecDot(dir12, dir12);
	if(fabs(a) < SMALL){
		//UG_LOG("a == 0: " << v1 << v2 << v3 << endl);
		return false;
	}

	number b = VecDot(dir12, proj);
	if(fabs(b) < SMALL){
		//UG_LOG("b == 0: " << v1 << v2 << v3 << endl);
		return false;
	}

	vector3 n1;
	VecScaleAdd(n1, -b / a, dir12, 1., proj);

	vector3 c1;
	VecScaleAdd(c1, 0.5, v1, 0.5, v2);

//	we also have to construct a plane. The planes normal is dir13.
	vector3 c2;
	VecScaleAdd(c2, 0.5, v1, 0.5, v3);

//	finally calculate the intersection of the line with the plane.
	number t;
	bool retVal = RayPlaneIntersection(centerOut, t, c1, n1, c2, dir13);
/*
	number dist1 = fabs(VecDistanceSq(centerOut, v1) - VecDistanceSq(centerOut, v2));
	number dist2 = fabs(VecDistanceSq(centerOut, v1) - VecDistanceSq(centerOut, v3));
	if( retVal &&
		(dist1 > SMALL || dist2 > SMALL))
	{
		UG_LOG("Center distance mismatch!!!\n");
		UG_LOG("v: " << v1 << v2 << v3 << endl);
		UG_LOG("|v2-v1| = " << dist1 <<", |v3-v1| = " << dist2 << endl);
	}*/
	return retVal;
}

bool FindNormal(vector3& normOut, const vector3& v)
{
//	normalize v to avoid problems with oversized vectors.
	vector3 n;
	VecNormalize(n, v);
	
//TODO: check whether 0.7 is a good threshold
	const number dotThreshold = 0.7;
	
//	try projections of the unit-normals onto the
//	plane which is defined by normOut and the origin.
	for(int i = 0; i < 3; ++i){
		vector3 e(0, 0, 0);
		e[i] = 1;
		number d = VecDot(e, n);
		if(fabs(d) < dotThreshold){
		//	the projection will be sufficient to calculate a normal.
			VecScale(n, n, d);
			VecSubtract(normOut, e, n);
			VecNormalize(normOut, normOut);
			return true;
		}
	}

	return false;
}

bool ConstructOrthonormalSystem(matrix33& matOut, const vector3& v,
								size_t vColInd)
{
	// the first expression has been always true for unsigned int (AV)
	if(/*vColInd < 0 ||*/ vColInd > 2)
		return false;
		
	vector3 newCols[2];
	vector3 n;
	
//	normalize v
	VecNormalize(n, v);
	
//	find a normal to n
	if(!FindNormal(newCols[0], n))
		return false;
		
//	calculate the last col
	VecCross(newCols[1], n, newCols[0]);
	
//	copy cols to matOut
	int nColCount = 0;
	for(size_t j = 0; j < 3; ++j){
		if(j == vColInd){
			for(size_t i = 0; i < 3; ++ i)
				matOut[i][j] = n[i];
		}
		else{
			for(size_t i = 0; i < 3; ++ i)
				matOut[i][j] = newCols[nColCount][i];
			++nColCount;
		}
	}
	
	return true;
}

void CalculateCovarianceMatrix(matrix33& matOut, const vector3* pointSet,
							  const vector3& center, size_t numPoints)
{
//	set all matrix entries to 0
	for(size_t i = 0; i < 3; ++i){
		for(size_t j = 0; j < 3; ++j)
			matOut[i][j] = 0;
	}
	
//	sum the vector-products
	for(size_t pInd = 0; pInd < numPoints; ++pInd){
		for(size_t i = 0; i < 3; ++i){
			for(size_t j = 0; j < 3; ++j){
				matOut[i][j] += (pointSet[pInd][i] - center[i]) * 
								 (pointSet[pInd][j] - center[j]);
			}
		}
	}
}

bool FindClosestPlane(vector3& centerOut, vector3& normalOut,
					  const vector3* pointSet, size_t numPoints)
{
//	calculate the center of the point set
	vector3 center;
	CalculateCenter(center, pointSet, numPoints);

//	calculate the covariance matrix of the point set
	matrix33 matCo;
	CalculateCovarianceMatrix(matCo, pointSet, center, numPoints);

//	find the eigenvector of smallest eigenvalue of the covariance matrix
	number lambda[3];
	vector3 ev[3];
	if(!CalculateEigenvalues(matCo, lambda[0], lambda[1], lambda[2],
							 ev[0], ev[1], ev[2]))
	{
		return false;
	}

	VecNormalize(normalOut, ev[0]);
	centerOut = center;

	return true;
}

bool TransformPointSetTo2D(vector2* pointSetOut, const vector3* pointSet,
						  size_t numPoints)
{
//	calculate the center of the point set
	vector3 center;
	vector3 normal;
	
	if(!FindClosestPlane(center, normal, pointSet, numPoints))
		return false;
	
//	lambda[0] is the smallest (absolute) eigenvalue.
//	ev[0] can be regarded as the normal of a plane through center.
//	we now have to find the matrix that transforms this plane
//	to a plane parallel to the x-y-plane.
	matrix33 matOrtho;
	if(!ConstructOrthonormalSystem(matOrtho, normal, 2))
		return false;
	
//	we need the inverse of this matrix. Since it is a orthonormal
//	matrix, the inverse corresponds to the transposed matrix.
//	move the point set and rotate it afterwards.
	for(size_t i = 0; i < numPoints; ++i){
		vector3 vTmpIn;
		vector3 vTmpOut;
		VecSubtract(vTmpIn, pointSet[i], center);
		TransposedMatVecMult(vTmpOut, matOrtho, vTmpIn);
		
	//	trivial projection
		pointSetOut[i].x() = vTmpOut.x();
		pointSetOut[i].y() = vTmpOut.y();
	}

//	done. return true.
	return true;
}

////////////////////////////////////////////////////////////////////////
bool RayRayIntersection3d(vector3& aOut, vector3& bOut,
						  const vector3& p0, const vector3& dir0,
						  const vector3& p1, const vector3& dir1)
{
	vector3 vNear, vAB;
	bool trueIntersection;
	vector3 q0, q1;
	VecAdd(q0, p0, dir0);
	VecAdd(q1, p1, dir1);

	IntersectLineSegments(p0.x(), p0.y(), p0.z(), q0.x(), q0.y(), q0.z(),
						  p1.x(), p1.y(), p1.z(), q1.x(), q1.y(), q1.z(),
						  true, SMALL,
						  aOut.x(), aOut.y(), aOut.z(), bOut.x(), bOut.y(), bOut.z(),
						  vNear.x(), vNear.y(), vNear.z(), vAB.x(), vAB.y(), vAB.z(),
						  trueIntersection);

	return trueIntersection;
}

////////////////////////////////////////////////////////////////////////
bool LineLineIntersection3d(vector3& aOut, vector3& bOut,
							const vector3& a1, const vector3& a2,
						  	const vector3& b1, const vector3& b2)
{
	vector3 vNear, vAB;
	bool trueIntersection;

	IntersectLineSegments(a1.x(), a1.y(), a1.z(), a2.x(), a2.y(), a2.z(),
						  b1.x(), b1.y(), b1.z(), b2.x(), b2.y(), b2.z(),
						  false, SMALL,
						  aOut.x(), aOut.y(), aOut.z(), bOut.x(), bOut.y(), bOut.z(),
						  vNear.x(), vNear.y(), vNear.z(), vAB.x(), vAB.y(), vAB.z(),
						  trueIntersection);

	return trueIntersection;
}

number DistanceLineToLine(const vector3& a1, const vector3& a2,
						  const vector3& b1, const vector3& b2)
{
	vector3 vA, vB, vNear, vAB;
	bool trueIntersection;

	IntersectLineSegments(a1.x(), a1.y(), a1.z(), a2.x(), a2.y(), a2.z(),
						  b1.x(), b1.y(), b1.z(), b2.x(), b2.y(), b2.z(),
						  false, SMALL,
						  vA.x(), vA.y(), vA.z(), vB.x(), vB.y(), vB.z(),
						  vNear.x(), vNear.y(), vNear.z(), vAB.x(), vAB.y(), vAB.z(),
						  trueIntersection);

	return VecLength(vAB);
}


bool RayCylinderIntersection(number& tMinOut, number& tMaxOut, const vector3& rayFrom,
							 const vector3& rayDir, const vector3& cylCenter,
							 const vector3& cylAxis, number cylRadius)
{
//	find the closest points on the ray and the cylinder axis
	vector3 vr, va;
	RayRayIntersection3d(vr, va, rayFrom, rayDir, cylCenter, cylAxis);

	number rayAxisDist = VecDistance(vr, va);

	if(rayAxisDist > cylRadius)
		return false;

//	parallel and normalized orthogonal component of rayDir regarding cylAxis
	vector3 dirParallel, dirOrtho;
	vector3 axisNormized;
	VecNormalize(axisNormized, cylAxis);
	VecScale(dirParallel, cylAxis, VecDot(axisNormized, rayDir));
	VecSubtract(dirOrtho, rayDir, dirParallel);

//	calculate the distance from rayFrom, where rayFrom+t*dirOrtho intersects the cylinder
	number orthoIntersectionDist = sqrt(cylRadius*cylRadius - rayAxisDist*rayAxisDist);

	number orthoLen = VecLength(dirOrtho);
	if(orthoLen < SMALL)
		return false;

//	the factor by which rayDir has to be scaled to reach the intersection point
//	starting from vr.
	number scaleFac = orthoIntersectionDist / orthoLen;

//	get the parameter at which vr lies on the given ray
	number t = 0;
	if(fabs(rayDir.x()) >= fabs(rayDir.y())){
		if(fabs(rayDir.x()) >= fabs(rayDir.z())){
			if(rayDir.x() == 0)
				return false;
			t = (vr.x() - rayFrom.x()) / rayDir.x();
		}
		else
			t = (vr.z() - rayFrom.z()) / rayDir.z();
	}
	else if(fabs(rayDir.y()) >= fabs(rayDir.z()))
		t = (vr.y() - rayFrom.y()) / rayDir.y();
	else
		t = (vr.z() - rayFrom.z()) / rayDir.z();

	tMinOut = t - scaleFac;
	tMaxOut = t + scaleFac;
	return true;
}


////////////////////////////////////////////////////////////////////////////////
number CalculateTetrahedronVolume(const vector3& a, const vector3& b,
								  const vector3& c, const vector3& d)
{
//
//	Assume a tetrahedron with vertices a, b, c, d, then the volume is given by
//
//	V = 1/6 * |VecDot( (a-d) , VecCross((b-d), (c-d)) )|
//
	number tetrahedronVolume;
	vector3 ad;
	vector3 bd;
	vector3 cd;
	vector3 cross;

	VecSubtract(ad, a, d);
	VecSubtract(bd, b, d);
	VecSubtract(cd, c, d);

	VecCross(cross, bd, cd);

	tetrahedronVolume = abs(VecDot(ad, cross) / 6);

	return tetrahedronVolume;
}

number CalculatePyramidVolume(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e)
{
	vector3 center, h_, h, c1, c2, da, ba, cb, cd;
	VecSubtract(da, d, a);
	VecSubtract(ba, b, a);
	VecSubtract(cb, c, b);
	VecSubtract(cd, c, d);
	VecCross(c1, da, ba);
	VecCross(c2, cb, cd);
	number A = 0.5 * VecLength(c1) + 0.5 * VecLength(c2);

	vector3 arr[] = { a, b, c, d };
	CalculateCenter(center, arr, 4);

	number height = DistancePointToPlane(e, center, c1);

	return 1.0 / 3.0 * A * height;
}

number CalculatePrismVolume(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e,
		const vector3& f) {
	number result = 0;
	vector3 center;
	vector3 arr[] = { a, b, c, d, e, f };
	CalculateCenter(center, arr, 6);

	result += CalculateTetrahedronVolume(a, b, c, center);
	result += CalculateTetrahedronVolume(d, e, f, center);

	result +=CalculatePyramidVolume(a, b, e, d, center);
	result +=CalculatePyramidVolume(b, c, f, e, center);
	result +=CalculatePyramidVolume(c, a, d, f, center);

	return result;
}

number CalculateHexahedronVolume(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e, const vector3& f,
		const vector3& g, const vector3& h) {
	number result = 0;
	vector3 arr[] = { a, b, c, d, e, f, g, h };
	vector3 center;
	CalculateCenter(center, arr, 8);

	// top and bottom
	result += CalculatePyramidVolume(a, b, c, d, center);
	result += CalculatePyramidVolume(e, f, g, h, center);

	// sides
	result += CalculatePyramidVolume(a, b, f, e, center);
	result += CalculatePyramidVolume(b, c, g, f, center);

	result += CalculatePyramidVolume(c, d, g, h, center);
	result += CalculatePyramidVolume(a, d, h, e, center);

	return result;
}

number CalculateOctahedronVolume(const vector3& a, const vector3& b,
		const vector3& c, const vector3& d, const vector3& e,
		const vector3& f) {
	number result = 0;

	MathVector<3> x31, x42, x51, n;
	MathVector<3> 			x01;

	VecSubtract(x31, d, b);
	VecSubtract(x42, e, c);
	VecSubtract(x51, f, b);
	VecCross(n, x31, x42);

	//VecSubtract(x31, d, b);
	//VecSubtract(x42, e, c);
	VecSubtract(x01, a, b);
	//VecCross(n, x31, x42);

	number volTopPyr 	= (1./6.) * fabs(VecDot(n, x51));
	number volBottomPyr = (1./6.) * fabs(VecDot(n, x01));

	result += volTopPyr + volBottomPyr;

	return result;
}

////////////////////////////////////////////////////////////////////////////////
//	Returns the BinomialCoefficient
int BinomCoeff(int n, int k)
{
//	if n == 0, we define result to be always zero iff k != 0
	if(!n&&k) return 0;

//	if denominator is greater than nominator: flip
	if(n-k>k) return BinomCoeff(n,n-k);

//	if equal binomCoeff is always one
	if(n==k) return 1;

//	do recursion
	return BinomCoeff(n-1,k)*n/(n-k);
}

}//	end of namespace

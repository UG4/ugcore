/**
 * \file math_util.cpp
 */

//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m05 d18

#include <vector>
#include "math_util.h"
#include "../ugmath.h"
#include "../ugmath_types.h"


using namespace std;

namespace ug
{

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
		if(d < dotThreshold){
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
	if(vColInd < 0 || vColInd > 2)
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
		pointSetOut[i].x = vTmpOut.x;
		pointSetOut[i].y = vTmpOut.y;
	}

//	done. return true.
	return true;
}

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

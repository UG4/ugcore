/*
 * finite_volume_util.h
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__
#define __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"

namespace ug{

///////////////////
// Volume of SCV
///////////////////

template <int TRefDim, int TWorldDim>
number VolumeOfSCV(const std::vector<MathVector<TWorldDim> >& vPoints);

/////////////////////////////////////////////
// Specialization for 1D Reference Dimension
template <>
number VolumeOfSCV<1,1>(const std::vector<MathVector<1> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 2, "Must be a line.");
	return VecDistance(vPoints[0], vPoints[1]);
}
template <>
number VolumeOfSCV<1,2>(const std::vector<MathVector<2> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 2, "Must be a line.");
	return VecDistance(vPoints[0], vPoints[1]);
}
template <>
number VolumeOfSCV<1,3>(const std::vector<MathVector<3> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 2, "Must be a line.");
	return VecDistance(vPoints[0], vPoints[1]);
}

/////////////////////////////////////////////
// Specialization for 2D Reference Dimension

template <>
number VolumeOfSCV<2,2>(const std::vector<MathVector<2> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 4, "Must be a quadrilateral.");

	const number tmp = (vPoints[3][1]-vPoints[1][1])*(vPoints[2][0]-vPoints[0][0])
				-(vPoints[3][0]-vPoints[1][0])*(vPoints[2][1]-vPoints[0][1]);
	return 0.5 * fabs( tmp );
}

template <>
number VolumeOfSCV<2,3>(const std::vector<MathVector<3> >& vPoints)
{
	UG_ASSERT(0, "Not implemented.");
	return -1;
}

/////////////////////////////////////////////
// Specialization for 3D Reference Dimension

template <>
number VolumeOfSCV<3,3>(const std::vector<MathVector<3> >& vPoints)
{
	UG_ASSERT(0, "Not implemented.");
	return -1;
}


///////////////////
// Normal on SCVF
///////////////////

template <int TRefDim, int TWorldDim>
void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const std::vector<MathVector<TWorldDim> >& vPoints)
{
	UG_ASSERT(0, "Not implemented.");
}

/////////////////////////////////////////////
// Specialization for 1D Reference Dimension

template <>
void NormalOnSCVF<1,1>(MathVector<1>& outNormal, const std::vector<MathVector<1> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 1, "Must be a point.");

	outNormal[0] = 1.0;
}

template <>
void NormalOnSCVF<1,2>(MathVector<2>& outNormal, const std::vector<MathVector<2> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 1, "Must be a point.");

	outNormal[0] = 1.0;
}

/////////////////////////////////////////////
// Specialization for 2D Reference Dimension

template <>
void NormalOnSCVF<2,2>(MathVector<2>& outNormal, const std::vector<MathVector<2> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 2, "Must be a line.");

	MathVector<2> diff = vPoints[1]; // center of element
	diff -= vPoints[0]; // edge midpoint

	outNormal[0] = diff[1];
	outNormal[1] = -diff[0];
}

/////////////////////////////////////////////
// Specialization for 3D Reference Dimension

template <>
void NormalOnSCVF<3,3>(MathVector<3>& outNormal, const std::vector<MathVector<3> >& vPoints)
{
	UG_ASSERT(vPoints.size() == 4, "Must be a quadrilateral.");

	MathVector<3> a, b;
	VecSubtract(a, vPoints[2], vPoints[0]);
	VecSubtract(b, vPoints[3], vPoints[1]);
	VecCross(outNormal, a,b);
	VecScale(outNormal, outNormal, 0.5);
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__ */

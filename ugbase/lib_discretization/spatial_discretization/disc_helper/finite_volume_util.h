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

/// averages positions by arithmetic mean
/** Arithmetic Mean of Positions
 * returns the arithmetic mean of positions
 *
 * \param[in]  vCornerCoords	positions
 * \param[in]  num				number of positions
 * \param[out] vOut				arithmetic mean of positions
 */
template <typename TPosition>
void AveragePositions(TPosition& vOut, const TPosition* vCornerCoords, size_t num)
{
	vOut = vCornerCoords[0];
	for(size_t j = 1; j < num; ++j)
	{
		vOut += vCornerCoords[j];
	}
	vOut *= 1./(number)num;
}



//////////////////////////
// Finite Volume Traits
//////////////////////////

template <typename TRefElem, int TWorldDim> struct finite_volume_traits;

/////////////////////////
// 1D Reference Element
/////////////////////////

template <> struct finite_volume_traits<ReferenceEdge, 1>
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;

	static void NormalOnSCVF(MathVector<1>& outNormal, const MathVector<1>* vCornerCoords)
		{ElementNormal<ReferenceVertex, 1>(outNormal, vCornerCoords);}

	typedef ReferenceEdge scv_type;
};

template <> struct finite_volume_traits<ReferenceEdge, 2>
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;

	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceEdge scv_type;
};

template <> struct finite_volume_traits<ReferenceEdge, 3>
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceEdge scv_type;
};

/////////////////////////
// 2D Reference Element
/////////////////////////

template <> struct finite_volume_traits<ReferenceTriangle, 2>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}

	typedef ReferenceQuadrilateral scv_type;
};

template <> struct finite_volume_traits<ReferenceTriangle, 3>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceQuadrilateral scv_type;
};

template <> struct finite_volume_traits<ReferenceQuadrilateral, 2>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}

	typedef ReferenceQuadrilateral scv_type;
};

template <> struct finite_volume_traits<ReferenceQuadrilateral, 3>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceQuadrilateral scv_type;
};

/////////////////////////
// 3D Reference Element
/////////////////////////

template <> struct finite_volume_traits<ReferenceTetrahedron, 3>
{
	const static size_t NumCornersOfSCVF = 4;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}

	typedef ReferenceHexahedron scv_type;
};

template <> struct finite_volume_traits<ReferencePrism, 3>
{
	const static size_t NumCornersOfSCVF = 4;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}

	typedef ReferenceHexahedron scv_type;
};

template <> struct finite_volume_traits<ReferencePyramid, 3>
{
	const static size_t NumCornersOfSCVF = 4;
	const static size_t MaxNumCornersOfSCV = 10;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}

	typedef ReferenceHexahedron scv_type;
};

template <> struct finite_volume_traits<ReferenceHexahedron, 3>
{
	const static size_t NumCornersOfSCVF = 4;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}

	typedef ReferenceHexahedron scv_type;
};


//////////////////////////
// Hanging Finite Volume Traits
//////////////////////////

template <typename TRefElem, int TWorldDim> struct hanging_finite_volume_traits;

/////////////////////////
// 1D Reference Element
/////////////////////////

template <> struct hanging_finite_volume_traits<ReferenceEdge, 1>
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;

	static void NormalOnSCVF(MathVector<1>& outNormal, const MathVector<1>* vCornerCoords)
		{ElementNormal<ReferenceVertex, 1>(outNormal, vCornerCoords);}

	typedef ReferenceEdge scv_type;
};

template <> struct hanging_finite_volume_traits<ReferenceEdge, 2>
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;

	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceEdge scv_type;
};

template <> struct hanging_finite_volume_traits<ReferenceEdge, 3>
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceEdge scv_type;
};

/////////////////////////
// 2D Reference Element
/////////////////////////

template <> struct hanging_finite_volume_traits<ReferenceTriangle, 2>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}

	typedef ReferenceQuadrilateral scv_type;
};

template <> struct hanging_finite_volume_traits<ReferenceTriangle, 3>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceQuadrilateral scv_type;
};

template <> struct hanging_finite_volume_traits<ReferenceQuadrilateral, 2>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}

	typedef ReferenceQuadrilateral scv_type;
};

template <> struct hanging_finite_volume_traits<ReferenceQuadrilateral, 3>
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_ASSERT(0, "Not implemented");}

	typedef ReferenceQuadrilateral scv_type;
};

/////////////////////////
// 3D Reference Element
/////////////////////////

template <> struct hanging_finite_volume_traits<ReferenceTetrahedron, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hanging_finite_volume_traits<ReferencePrism, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hanging_finite_volume_traits<ReferencePyramid, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 10;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hanging_finite_volume_traits<ReferenceHexahedron, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};


///////////////////
// Functions
///////////////////

template <typename TRefElem, int TWorldDim>
void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords)
{
	finite_volume_traits<TRefElem, TWorldDim>::NormalOnSCVF(outNormal, vCornerCoords);
}

template <typename TRefElem, int TWorldDim>
void HangingNormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords)
{
	hanging_finite_volume_traits<TRefElem, TWorldDim>::NormalOnSCVF(outNormal, vCornerCoords);
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPACIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__ */

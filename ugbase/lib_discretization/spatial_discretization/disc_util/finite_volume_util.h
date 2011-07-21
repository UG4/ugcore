/*
 * finite_volume_util.h
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__
#define __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"

#include "lib_discretization/common/geometry_util.h"

namespace ug{

/// averages positions by arithmetic mean
/**
 * Arithmetic Mean of Positions
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



////////////////////////////////////////////////////////////////////////////////
// Finite Volume Traits
////////////////////////////////////////////////////////////////////////////////

/// Traits for Finite Volumes (dummy implementation)
template <typename TRefElem, int TWorldDim> struct fv1_traits
{
	const static size_t NumCornersOfSCVF;
	const static size_t MaxNumCornersOfSCV;

	static void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords);

	typedef void scv_type;
};

/////////////////////////
// 1D Reference Element
/////////////////////////

// \todo: NumCornersOfSCVF = 2 is to much, handle after checking
struct fv1_traits_ReferenceEdge
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 2;
	typedef ReferenceEdge scv_type;
};

template <> struct fv1_traits<ReferenceEdge, 1> : public fv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<1>& outNormal, const MathVector<1>* vCornerCoords)
		{ElementNormal<ReferenceVertex, 1>(outNormal, vCornerCoords);}
};

template <> struct fv1_traits<ReferenceEdge, 2> : public fv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{VecSubtract(outNormal, vCornerCoords[1], vCornerCoords[0]);}
};

template <> struct fv1_traits<ReferenceEdge, 3> : public fv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{throw(UGFatalError("Not implemented"));}
};

/////////////////////////
// 2D Reference Element
/////////////////////////

struct fv1_traits_ReferenceFace
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;
	typedef ReferenceQuadrilateral scv_type;
};

struct fv1_traits_ReferenceFace2d : public fv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}
};

struct fv1_traits_ReferenceFace3d : public fv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{throw(UGFatalError("Not implemented"));}
};

template <> struct fv1_traits<ReferenceTriangle, 2> : public fv1_traits_ReferenceFace2d{};
template <> struct fv1_traits<ReferenceTriangle, 3> : public fv1_traits_ReferenceFace3d{};

template <> struct fv1_traits<ReferenceQuadrilateral, 2> : public fv1_traits_ReferenceFace2d{};
template <> struct fv1_traits<ReferenceQuadrilateral, 3> : public fv1_traits_ReferenceFace3d{};

/////////////////////////
// 3D Reference Element
/////////////////////////

struct fv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 4;
	typedef ReferenceHexahedron scv_type;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}
};

template <> struct fv1_traits<ReferenceTetrahedron, 3> : public fv1_traits_ReferenceVolume
{
	const static size_t MaxNumCornersOfSCV = 8;

};

template <> struct fv1_traits<ReferencePrism, 3> : public fv1_traits_ReferenceVolume
{
	const static size_t MaxNumCornersOfSCV = 8;
};

template <> struct fv1_traits<ReferencePyramid, 3> : public fv1_traits_ReferenceVolume
{
	const static size_t MaxNumCornersOfSCV = 10;
};

template <> struct fv1_traits<ReferenceHexahedron, 3> : public fv1_traits_ReferenceVolume
{
	const static size_t MaxNumCornersOfSCV = 8;
};


////////////////////////////////////////////////////////////////////////////////
// Hanging Finite Volume Traits
////////////////////////////////////////////////////////////////////////////////

///	Traits for hanging finite volume (dummy implementation)
template <typename TRefElem, int TWorldDim> struct hfv1_traits
{
	const static size_t NumCornersOfSCVF;
	const static size_t MaxNumCornersOfSCV;

	static void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords);

	typedef void scv_type;
};

/////////////////////////
// 1D Reference Element
/////////////////////////

struct hfv1_traits_ReferenceEdge
{
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;
	typedef ReferenceEdge scv_type;
};

template <> struct hfv1_traits<ReferenceEdge, 1> : public hfv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<1>& outNormal, const MathVector<1>* vCornerCoords)
		{ElementNormal<ReferenceVertex, 1>(outNormal, vCornerCoords);}
};

template <> struct hfv1_traits<ReferenceEdge, 2> : public hfv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{throw(UGFatalError("Not implemented"));}
};

template <> struct hfv1_traits<ReferenceEdge, 3> : public hfv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{throw(UGFatalError("Not implemented"));}
};

/////////////////////////
// 2D Reference Element
/////////////////////////

struct hfv1_traits_ReferenceFace
{
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;
	typedef ReferenceQuadrilateral scv_type;
};

template <> struct hfv1_traits<ReferenceTriangle, 2> : public hfv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}
};

template <> struct hfv1_traits<ReferenceTriangle, 3> : public hfv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{throw(UGFatalError("Not implemented"));}
};

template <> struct hfv1_traits<ReferenceQuadrilateral, 2> : public hfv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}
};

template <> struct hfv1_traits<ReferenceQuadrilateral, 3> : public hfv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{throw(UGFatalError("Not implemented"));}
};

/////////////////////////
// 3D Reference Element
/////////////////////////

template <> struct hfv1_traits<ReferenceTetrahedron, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hfv1_traits<ReferencePrism, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hfv1_traits<ReferencePyramid, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 10;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hfv1_traits<ReferenceHexahedron, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};


////////////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////////////

template <typename TRefElem, int TWorldDim>
void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords)
{
	fv1_traits<TRefElem, TWorldDim>::NormalOnSCVF(outNormal, vCornerCoords);
}

template <typename TRefElem, int TWorldDim>
void HangingNormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords)
{
	hfv1_traits<TRefElem, TWorldDim>::NormalOnSCVF(outNormal, vCornerCoords);
}

////////////////////////////////////////////////////////////////////////////////
// Dimension dependent traits
////////////////////////////////////////////////////////////////////////////////

///	Traits for Finite Volumes in a dimension
template <int TDim, int TWorldDim> struct fv1_dim_traits
{
//	inherited from fve_traits
	const static size_t NumCornersOfSCVF;
	const static size_t MaxNumCornersOfSCV;
	static void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords);
	typedef void scv_type;

//	defined only for dimension
	static const size_t maxNumSCVF;
	static const size_t maxNumSCV;
	static const size_t maxNSH;

};

/////////////////////////
// 1D
/////////////////////////

struct fv1_dim_traits_ReferenceEdge
{
	static const size_t maxNumSCVF = 2;
	static const size_t maxNumSCV = 2;
	static const size_t maxNSH = maxNumSCV;
};

template <> struct fv1_dim_traits<1, 1>
	: public fv1_traits<ReferenceEdge, 1>, public fv1_dim_traits_ReferenceEdge {};

template <> struct fv1_dim_traits<1, 2>
	: public fv1_traits<ReferenceEdge, 2>, public fv1_dim_traits_ReferenceEdge {};

template <> struct fv1_dim_traits<1, 3>
 	: public fv1_traits<ReferenceEdge, 3>, public fv1_dim_traits_ReferenceEdge {};

/////////////////////////
// 2D
/////////////////////////

struct fv1_dim_traits_ReferenceFace
{
	static const size_t maxNumSCVF = 4;
	static const size_t maxNumSCV = 4;
	static const size_t maxNSH = maxNumSCV;
};

template <> struct fv1_dim_traits<2, 2>
	: public fv1_traits_ReferenceFace2d, public fv1_dim_traits_ReferenceFace {};

template <> struct fv1_dim_traits<2, 3>
	: public fv1_traits_ReferenceFace3d, public fv1_dim_traits_ReferenceFace {};

/////////////////////////
// 3D
/////////////////////////

struct fv1_dim_traits_ReferenceVolume
{
	static const size_t maxNumSCVF = 12;
	static const size_t maxNumSCV = 8;
	static const size_t maxNSH = maxNumSCV;
};


template <> struct fv1_dim_traits<3, 3>
	: public fv1_traits_ReferenceVolume, public fv1_dim_traits_ReferenceVolume
{
	const static size_t MaxNumCornersOfSCV = 10;
};


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__DISC_HELPER__FINITE_VOLUME_UTIL__ */

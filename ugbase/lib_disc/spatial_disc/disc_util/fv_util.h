
#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_UTIL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_UTIL__

// extern libraries
#include <cmath>
#include <vector>

// other ug4 modules
#include "common/common.h"

#include "lib_disc/common/geometry_util.h"
#include "lib_disc/local_finite_element/lagrange/lagrange.h"

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
//	maximum for dimension
	static const size_t maxNumSCVF;
	static const size_t maxNumSCV;
	static const size_t maxNSH;

//	number of corners of scvf
	const static size_t NumCornersOfSCVF;

//	maximum of corners of scv
	const static size_t NumCornersOfSCV;
	
//	maximum of corners of bf
	const static size_t NumCornersOfBF;

//	computes the normal to a scvf
	static void NormalOnSCVF(MathVector<TWorldDim>& outNormal,
	                         const MathVector<TWorldDim>* vSCVFCorner,
	                         const MathVector<TWorldDim>* vElemCorner);

//	computes the normal to a bf
	static void NormalOnBF(MathVector<TWorldDim>& outNormal,
							 const MathVector<TWorldDim>* vSCVFCorner,
							 const MathVector<TWorldDim>* vElemCorner);

//	types of scv and scvf and bf
	typedef void scv_type;
	typedef void scvf_type;
	typedef void bf_type;
};

/////////////////////////
// 1D Reference Element
/////////////////////////

struct fv1_traits_ReferenceEdge
{
	static const size_t maxNumSCVF = 1;
	static const size_t maxNumSCV = 2;
	static const size_t maxNSH = maxNumSCV;

	const static size_t NumCornersOfSCVF = 1;
	const static size_t NumCornersOfSCV = 2;
	const static size_t NumCornersOfBF = 1;

	typedef ReferenceEdge scv_type;
	typedef ReferenceVertex scvf_type;
	typedef ReferenceVertex bf_type;
};

template <> struct fv1_traits<ReferenceEdge, 1> : public fv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<1>& outNormal,
	                         const MathVector<1>* vSCVFCorner,
	                         const MathVector<1>* vElemCorner)
	{
		ElementNormal<ReferenceVertex, 1>(outNormal, vSCVFCorner);
		VecNormalize(outNormal, outNormal);
	}
	static void NormalOnBF(MathVector<1>& outNormal,
	                       const MathVector<1>* vSCVFCorner,
	                       const MathVector<1>* vElemCorner)
	{
		NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
	}
};

template <> struct fv1_traits<ReferenceEdge, 2> : public fv1_traits_ReferenceEdge
{
		static void NormalOnSCVF(MathVector<2>& outNormal,
		                         const MathVector<2>* vSCVFCorner,
		                         const MathVector<2>* vElemCorner)
		{
			VecSubtract(outNormal, vElemCorner[1], vElemCorner[0]);
			VecNormalize(outNormal, outNormal);
		}
		static void NormalOnBF(MathVector<2>& outNormal,
		                       const MathVector<2>* vSCVFCorner,
		                       const MathVector<2>* vElemCorner)
		{
			NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
		}
};

template <> struct fv1_traits<ReferenceEdge, 3> : public fv1_traits_ReferenceEdge
{
		static void NormalOnSCVF(MathVector<3>& outNormal,
		                         const MathVector<3>* vSCVFCorner,
		                         const MathVector<3>* vElemCorner)
		{
			VecSubtract(outNormal, vElemCorner[1], vElemCorner[0]);
			VecNormalize(outNormal, outNormal);
		}
		static void NormalOnBF(MathVector<3>& outNormal,
		                       const MathVector<3>* vSCVFCorner,
		                       const MathVector<3>* vElemCorner)
		{
			NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
		}
};

/////////////////////////
// 2D Reference Element
/////////////////////////

struct fv1_traits_ReferenceFace
{
	static const size_t maxNumSCVF = 4;
	static const size_t maxNumSCV = 4;
	static const size_t maxNSH = maxNumSCV;

	const static size_t NumCornersOfSCVF = 2;
	const static size_t NumCornersOfSCV = 4;
	const static size_t NumCornersOfBF = 2;
	
	typedef ReferenceQuadrilateral scv_type;
	typedef ReferenceEdge scvf_type;
	typedef ReferenceEdge bf_type;
};

struct fv1_traits_ReferenceFace2d : public fv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<2>& outNormal,
							 const MathVector<2>* vSCVFCorner,
							 const MathVector<2>* vElemCorner)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vSCVFCorner);}
	static void NormalOnBF(MathVector<2>& outNormal,
	                       const MathVector<2>* vSCVFCorner,
	                       const MathVector<2>* vElemCorner)
	{
		NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
	}
};

struct fv1_traits_ReferenceFace3d : public fv1_traits_ReferenceFace
{
	template <typename TElem>
	static void NormalOnSCVF_Face(MathVector<3>& outNormal,
							 const MathVector<3>* vSCVFCorner,
							 const MathVector<3>* vElemCorner)
	{
		// compute Normal to Face (right-handed)
		MathVector<3> ElemNormal;
		ElementNormal<TElem,3>(ElemNormal, vElemCorner);

		// compute a Point such that the triangle given by vSCVFCorner and p
		// has the following property:
		// a) The new Triangle is normal to the ElementFace
		// b) The new Triangle contains the scvf
		MathVector<3> p;
		VecAdd(p, vSCVFCorner[0], ElemNormal);

		MathVector<3> vNewTriangleCorner[3];
		vNewTriangleCorner[0] = vSCVFCorner[0];
		vNewTriangleCorner[1] = vSCVFCorner[1];
		vNewTriangleCorner[2] = p;

		// now compute Normal to new triangle. This Normal will be normal to
		// the scvf and within the original element
		ElementNormal<ReferenceTriangle, 3>(outNormal, vNewTriangleCorner);

		// scale to normal to the size of the scvf
		const number size = VecDistance(vSCVFCorner[0], vSCVFCorner[1]);
		VecNormalize(outNormal, outNormal);
		VecScale(outNormal, outNormal, size);
	}
	static void NormalOnSCVF(MathVector<3>& outNormal,
							 const MathVector<3>* vSCVFCorner,
							 const MathVector<3>* vElemCorner)
		{
			UG_THROW("not implemented.")
		}
	static void NormalOnBF(MathVector<3>& outNormal,
	                       const MathVector<3>* vSCVFCorner,
	                       const MathVector<3>* vElemCorner)
		{
			UG_THROW("not implemented.")
		}
};

template <> struct fv1_traits<ReferenceTriangle, 2> : public fv1_traits_ReferenceFace2d{};
template <> struct fv1_traits<ReferenceTriangle, 3> : public fv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<3>& outNormal,
							 const MathVector<3>* vSCVFCorner,
							 const MathVector<3>* vElemCorner)
		{fv1_traits_ReferenceFace3d::NormalOnSCVF_Face<ReferenceTriangle>(outNormal, vSCVFCorner, vElemCorner);}
	static void NormalOnBF(MathVector<3>& outNormal,
	                       const MathVector<3>* vSCVFCorner,
	                       const MathVector<3>* vElemCorner)
	{
		NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
	}
};

template <> struct fv1_traits<ReferenceQuadrilateral, 2> : public fv1_traits_ReferenceFace2d{};
template <> struct fv1_traits<ReferenceQuadrilateral, 3> : public fv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<3>& outNormal,
							 const MathVector<3>* vSCVFCorner,
							 const MathVector<3>* vElemCorner)
		{fv1_traits_ReferenceFace3d::NormalOnSCVF_Face<ReferenceQuadrilateral>(outNormal, vSCVFCorner, vElemCorner);}
	static void NormalOnBF(MathVector<3>& outNormal,
	                       const MathVector<3>* vSCVFCorner,
	                       const MathVector<3>* vElemCorner)
	{
		NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
	}
};

/////////////////////////
// 3D Reference Element
/////////////////////////

struct fv1_traits_ReferenceVolume
{
	static const size_t maxNumSCVF = 24;
	static const size_t maxNumSCV = 32;
	static const size_t maxNSH = 8;

	const static size_t NumCornersOfSCVF = 4;
	const static size_t NumCornersOfSCV = 8;
	const static size_t NumCornersOfBF = 4;

	typedef ReferenceHexahedron scv_type;
	typedef ReferenceQuadrilateral scvf_type;
	typedef ReferenceQuadrilateral bf_type;

	static void NormalOnSCVF(MathVector<3>& outNormal,
	                         const MathVector<3>* vSCVFCorner,
	                         const MathVector<3>* vElemCorner)
		{ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vSCVFCorner);}
	static void NormalOnBF(MathVector<3>& outNormal,
	                       const MathVector<3>* vSCVFCorner,
	                       const MathVector<3>* vElemCorner)
	{
		NormalOnSCVF(outNormal, vSCVFCorner, vElemCorner);
	}
};

template <> struct fv1_traits<ReferenceTetrahedron, 3> : public fv1_traits_ReferenceVolume{};
template <> struct fv1_traits<ReferencePrism, 3> : public fv1_traits_ReferenceVolume{};
template <> struct fv1_traits<ReferenceHexahedron, 3> : public fv1_traits_ReferenceVolume{};
template <> struct fv1_traits<ReferenceOctahedron, 3> : public fv1_traits_ReferenceVolume{};

// For Pyramids we use triangular scvf, since the quadrilateral scvf would not be
// flat in general by the positions where its corners are placed
template <> struct fv1_traits<ReferencePyramid, 3> : public fv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 3; // triangles
	const static size_t NumCornersOfSCV = 4;  // tetrahedrons
	const static size_t NumCornersOfBF = 4;   // quadrilaterals

	typedef ReferenceTetrahedron scv_type;
	typedef ReferenceTriangle scvf_type;
	typedef ReferenceQuadrilateral bf_type;

	static void NormalOnSCVF(MathVector<3>& outNormal,
							 const MathVector<3>* vSCVFCorner,
							 const MathVector<3>* vElemCorner)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vSCVFCorner);}
	static void NormalOnBF(MathVector<3>& outNormal,
	                       const MathVector<3>* vSCVFCorner,
	                       const MathVector<3>* vElemCorner)
	{
		ElementNormal<ReferenceQuadrilateral, 3>(outNormal, vSCVFCorner);
	}
};

////////////////////////////////////////////////////////////////////////////////
// Dimension dependent traits DIM FV1
////////////////////////////////////////////////////////////////////////////////

///	Traits for Finite Volumes in a dimension
template <int TDim, int TWorldDim> struct fv1_dim_traits;

template <> struct fv1_dim_traits<1, 1> : public fv1_traits<ReferenceEdge, 1> {};
template <> struct fv1_dim_traits<1, 2> : public fv1_traits<ReferenceEdge, 2> {};
template <> struct fv1_dim_traits<1, 3> : public fv1_traits<ReferenceEdge, 3> {};

template <> struct fv1_dim_traits<2, 2> : public fv1_traits_ReferenceFace2d {};
template <> struct fv1_dim_traits<2, 3> : public fv1_traits_ReferenceFace3d {
	static void NormalOnSCVF(MathVector<3>& outNormal,
							 const MathVector<3>* vSCVFCorner,
							 const MathVector<3>* vElemCorner)
	{
		// Little bit dirty, but should be correct:
		// Even if the true element has more than three vertices (quadrilateral),
		// we only need three to compute the direction of the normal ElemNormal in NormalOnSCVF_Face,
		// the norm is not needed!
		fv1_traits_ReferenceFace3d::NormalOnSCVF_Face<ReferenceTriangle>(outNormal, vSCVFCorner, vElemCorner);
		//UG_THROW("Not implemented.")
	}

};

template <> struct fv1_dim_traits<3, 3>	: public fv1_traits_ReferenceVolume {};

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
		{UG_THROW("Not implemented");}
};

template <> struct hfv1_traits<ReferenceEdge, 3> : public hfv1_traits_ReferenceEdge
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_THROW("Not implemented");}
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
		{UG_THROW("Not implemented");}
};

template <> struct hfv1_traits<ReferenceQuadrilateral, 2> : public hfv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<2>& outNormal, const MathVector<2>* vCornerCoords)
		{ElementNormal<ReferenceEdge, 2>(outNormal, vCornerCoords);}
};

template <> struct hfv1_traits<ReferenceQuadrilateral, 3> : public hfv1_traits_ReferenceFace
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{UG_THROW("Not implemented");}
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
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hfv1_traits<ReferenceOctahedron, 3>
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;

	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
		{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <int TDim> struct hdimfv1_traits
{
	typedef void scv_type;
	typedef void elem_type_0;
	typedef void elem_type_1;
	typedef void elem_type_2;
	typedef void elem_type_3;
	typedef void elem_type_4;
	const static size_t NumCornersOfSCVF;
	const static size_t MaxNumCornersOfSCV;
};

template <> struct hdimfv1_traits<1>
{
	typedef ReferenceEdge scv_type;
	typedef ReferenceEdge elem_type_0;
	typedef ReferenceEdge elem_type_1;
	typedef ReferenceEdge elem_type_2;
	typedef ReferenceEdge elem_type_3;
	typedef ReferenceEdge elem_type_4;
	const static size_t NumCornersOfSCVF = 1;
	const static size_t MaxNumCornersOfSCV = 2;
};

template <> struct hdimfv1_traits<2>
{
	typedef ReferenceQuadrilateral scv_type;
	typedef ReferenceTriangle elem_type_0;
	typedef ReferenceQuadrilateral elem_type_1;
	typedef ReferenceTriangle elem_type_2;
	typedef ReferenceQuadrilateral elem_type_3;
	typedef ReferenceQuadrilateral elem_type_4;
	const static size_t NumCornersOfSCVF = 2;
	const static size_t MaxNumCornersOfSCV = 4;
};

template <> struct hdimfv1_traits<3>
{
	typedef ReferenceTetrahedron scv_type;
	typedef ReferenceTetrahedron elem_type_0;
	typedef ReferencePyramid elem_type_1;
	typedef ReferencePrism elem_type_2;
	typedef ReferenceHexahedron elem_type_3;
	typedef ReferenceOctahedron elem_type_4;
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 10;
};


////////////////////////////////////////////////////////////////////////////////
// Functions
////////////////////////////////////////////////////////////////////////////////

template <typename TRefElem, int TWorldDim>
void HangingNormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords)
{
	hfv1_traits<TRefElem, TWorldDim>::NormalOnSCVF(outNormal, vCornerCoords);
}

////////////////////////////////////////////////////////////////////////////////
// FVHO: Finite Volume for Higher Order traits
////////////////////////////////////////////////////////////////////////////////

///	Traits for Finite Volumes of higher order
template <int TOrder, typename TRefElem, int TWorldDim> struct fvho_traits
{
//	can be inherited from fv1_traits (since the same)
	const static size_t NumCornersOfSCVF;
	const static size_t MaxNumCornersOfSCV;
	static void NormalOnSCVF(MathVector<TWorldDim>& outNormal, const MathVector<TWorldDim>* vCornerCoords);
	typedef void scv_type;

//	own data
	const static size_t NumSubElem;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FINITE_VOLUME_UTIL__ */

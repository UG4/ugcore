/*
 * Copyright (c) 2010-2021:  G-CSC, Goethe University Frankfurt
 * Authors: Andreas Vogel, Dmitrij Logashenko, Martin Stepniewski
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

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

///	Base class, some fields are redefined in the instantiations for particular elements
template <typename TRefElem> struct fv1_traits_most_common
{
///	type of reference element
	typedef TRefElem ref_elem_type;

///	number of SubControlVolumes (for most of the elements - overridden for e.g. ROID_PYRAMID and ROID_OCTAHEDRON)
	static const size_t numSCV = ref_elem_type::numCorners;
	
///	number of SubControlVolumeFaces (for most of the elements - overridden for e.g. ROID_PYRAMID and ROID_OCTAHEDRON)
	static const size_t numSCVF = ref_elem_type::numEdges;

/// returns the 'from' and 'to' corner indices for a scvf
	static size_t scvf_from_to
	(
		const ref_elem_type& refElem, ///< the reference element object
		size_t i, ///< index of the scvf
		size_t ft ///< 0 = from, 1 = to
	)
	{
		return refElem.id(1, i, 0, ft);
	}
	
///	returns the node id for a scv
	static size_t scv_node_id
	(
		const ref_elem_type& refElem, ///< the reference element object
		size_t i ///< index of the scv
	)
	{
		return i;
	}
};

/// Traits for Finite Volumes (dummy implementation, s. the instantiations below)
template <typename TRefElem, int TWorldDim> struct fv1_traits
:	public fv1_traits_most_common<TRefElem>
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

template <> struct fv1_traits<ReferenceEdge, 1>
:	public fv1_traits_ReferenceEdge,
	public fv1_traits_most_common<ReferenceEdge>
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

template <> struct fv1_traits<ReferenceEdge, 2>
:	public fv1_traits_ReferenceEdge,
	public fv1_traits_most_common<ReferenceEdge>
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

template <> struct fv1_traits<ReferenceEdge, 3>
:	public fv1_traits_ReferenceEdge,
	public fv1_traits_most_common<ReferenceEdge>
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

template <> struct fv1_traits<ReferenceTriangle, 2>
:	public fv1_traits_ReferenceFace2d,
	public fv1_traits_most_common<ReferenceTriangle>
{};
template <> struct fv1_traits<ReferenceTriangle, 3>
:	public fv1_traits_ReferenceFace,
	public fv1_traits_most_common<ReferenceTriangle>
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

template <> struct fv1_traits<ReferenceQuadrilateral, 2>
:	public fv1_traits_ReferenceFace2d,
	public fv1_traits_most_common<ReferenceQuadrilateral>
{};
template <> struct fv1_traits<ReferenceQuadrilateral, 3>
:	public fv1_traits_ReferenceFace,
	public fv1_traits_most_common<ReferenceQuadrilateral>
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

template <> struct fv1_traits<ReferenceTetrahedron, 3>
:	public fv1_traits_ReferenceVolume,
	public fv1_traits_most_common<ReferenceTetrahedron>
{};
template <> struct fv1_traits<ReferencePrism, 3>
:	public fv1_traits_ReferenceVolume,
	public fv1_traits_most_common<ReferencePrism>
{};
template <> struct fv1_traits<ReferenceHexahedron, 3>
:	public fv1_traits_ReferenceVolume,
	public fv1_traits_most_common<ReferenceHexahedron>
{};

/// Pyramids: dimension-independent part of the FV1 traits
/**
 * For Pyramids we use triangular scvf, since the quadrilateral scvf would not be
 * flat in general by the positions where its corners are placed
 */
struct fv1_traits_ReferencePyramid
:	public fv1_traits_ReferenceVolume
	// Remark: Pyramid is a special case, fv1_traits_most_common<ReferencePyramid> is not inherited here!
{
	static const size_t numSCV = 4 * ReferencePyramid::numEdges; ///< overridden field from fv1_traits_most_common
	
	static const size_t numSCVF = 2 * ReferencePyramid::numEdges; ///< overridden field from fv1_traits_most_common

	const static size_t NumCornersOfSCVF = 3; // triangles
	const static size_t NumCornersOfSCV = 4;  // tetrahedrons
	const static size_t NumCornersOfBF = 4;   // quadrilaterals

	typedef ReferenceTetrahedron scv_type;
	typedef ReferenceTriangle scvf_type;
	typedef ReferenceQuadrilateral bf_type;

/// returns the 'from' and 'to' corner indices for a scvf (overridden function from fv1_traits_most_common)
	static size_t scvf_from_to
	(
		const ReferencePyramid& refElem, ///< the reference element object
		size_t i, ///< index of the scvf
		size_t ft ///< 0 = from, 1 = to
	)
	{ // map according to the order defined in ComputeSCVFMidID
		return refElem.id(1, i/2, 0, ft);
	}
	
///	returns the node id for a scv (overridden function from fv1_traits_most_common)
	static size_t scv_node_id
	(
		const ReferencePyramid& refElem, ///< the reference element object
		size_t i ///< index of the scv
	)
	{ // map according to order defined in ComputeSCVMidID
		if(i%2 == 0){
			return refElem.id(1, i/4, 0, 0); // from
		} else {
			return refElem.id(1, i/4, 0, 1); // to
		}
	}
};
/// Pyramids: the FV1 traits
/**
 * For Pyramids we use triangular scvf, since the quadrilateral scvf would not be
 * flat in general by the positions where its corners are placed
 */
template <> struct fv1_traits<ReferencePyramid, 3> : public fv1_traits_ReferencePyramid
{
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

/// Octahedra: dimension-independent part of the FV1 traits
struct fv1_traits_ReferenceOctahedron
:	public fv1_traits_ReferenceVolume
	// Remark: Octahedron is a special case, fv1_traits_most_common<ReferenceOctahedron> is not inherited here!
{
/**
 * The octahedral reference element contains an implicit interior
 * substructure that is constructed by several geometric objects,
 * i.e. imaginary subfaces (8 triangles), subedges (2) and subvolumes (4 tetrahedra)
 * resulting from the division into 4 subtetrahedra alongside inner edge 3->1.
 * The dual fv1 geometry consists of the original hexahedral SCVs in each of the
 * 4 subtetrahedra.
 */
	static const size_t numSCV = 16; ///< overridden field from fv1_traits_most_common
	
	static const size_t numSCVF = 24; ///< overridden field from fv1_traits_most_common
	//	Remark: Special case for octahedron, scvf not mappable by edges.
	
//	maximum dimension of substructure objects
	enum{MAXDIM = 3};

//	maximum number of substructure objects
	enum{MAXSUBSTRUCTOBJECTS = 8};

//	maximum number of substructure corners
	enum{MAXSUBSTRUCTCORNERS = 4};

/// returns the id of corner j of obj i in dimension dim_i of the implicit interior substructure
/**
 * The octahedral reference element contains an implicit interior
 * substructure that is constructed by several geometric objects, that
 * are mapped by a reference element by themselves. This method returns the
 * id (w.r.t. this reference element) of a sub-geometric object that is
 * part of a sub-geometric object of the implicit interior substructure of
 * this reference element.
 *
 * \param[in]	dim_i		dimension of sub geometric object
 * \param[in]	i			id of sub geometric object
 * \param[in]	j			number of obj contained in the sub-object
 * \returns		id of the j'th corner that is
 * 				contained in the i*th (sub-)geom object of dimension dim_i
 */
	static int substruct_coID(int dim_i, size_t i, size_t j)
	{
	// 	Corner indices of implicit interior substructure Geometric Objects
		static int substruct_coID[MAXDIM+1][MAXSUBSTRUCTOBJECTS][MAXSUBSTRUCTCORNERS];

		// subedge 0 = (3,1)
		substruct_coID[EDGE][0][0] = 3;
		substruct_coID[EDGE][0][1] = 1;
		// subedge 1 = (1,3)
		substruct_coID[EDGE][1][0] = 1;
		substruct_coID[EDGE][1][1] = 3;

		// subface 0 = (1,2,3)
		substruct_coID[FACE][0][0] = 1;
		substruct_coID[FACE][0][1] = 2;
		substruct_coID[FACE][0][2] = 3;
		// subface 1 = (1,3,2)
		substruct_coID[FACE][1][0] = 1;
		substruct_coID[FACE][1][1] = 3;
		substruct_coID[FACE][1][2] = 2;
		// subface 2 = (1,3,4)
		substruct_coID[FACE][2][0] = 1;
		substruct_coID[FACE][2][1] = 3;
		substruct_coID[FACE][2][2] = 4;
		// subface 3 = (1,4,3)
		substruct_coID[FACE][3][0] = 1;
		substruct_coID[FACE][3][1] = 4;
		substruct_coID[FACE][3][2] = 3;
		// subface 4 = (1,3,5)
		substruct_coID[FACE][4][0] = 1;
		substruct_coID[FACE][4][1] = 3;
		substruct_coID[FACE][4][2] = 5;
		// subface 5 = (1,5,3)
		substruct_coID[FACE][5][0] = 1;
		substruct_coID[FACE][5][1] = 5;
		substruct_coID[FACE][5][2] = 3;
		// subface 6 = (1,0,3)
		substruct_coID[FACE][6][0] = 1;
		substruct_coID[FACE][6][1] = 0;
		substruct_coID[FACE][6][2] = 3;
		// subface 7 = (1,0,3)
		substruct_coID[FACE][7][0] = 1;
		substruct_coID[FACE][7][1] = 3;
		substruct_coID[FACE][7][2] = 0;

		// subvolume 0 = (1,2,3,5)
		substruct_coID[VOLUME][0][0] = 1;
		substruct_coID[VOLUME][0][1] = 2;
		substruct_coID[VOLUME][0][2] = 3;
		substruct_coID[VOLUME][0][3] = 5;
		// subvolume 1 = (1,3,4,5)
		substruct_coID[VOLUME][1][0] = 1;
		substruct_coID[VOLUME][1][1] = 3;
		substruct_coID[VOLUME][1][2] = 4;
		substruct_coID[VOLUME][1][3] = 5;
		// subvolume 2 = (1,2,3,0)
		substruct_coID[VOLUME][2][0] = 1;
		substruct_coID[VOLUME][2][1] = 2;
		substruct_coID[VOLUME][2][2] = 3;
		substruct_coID[VOLUME][2][3] = 0;
		// subvolume 3 = (1,3,4,0)
		substruct_coID[VOLUME][3][0] = 1;
		substruct_coID[VOLUME][3][1] = 3;
		substruct_coID[VOLUME][3][2] = 4;
		substruct_coID[VOLUME][3][3] = 0;

		return substruct_coID[dim_i][i][j];
	}

/// returns the number of implicit interior substructure geometric objects of dim
/**
 * The octahedral reference element contains an implicit interior
 * substructure that is constructed by several geometric objects, that
 * are mapped by a reference element by themselves. This method returns how
 * many (sub-)geometric objects of a given dimension are contained in the
 * implicit interior substructure of this reference element.
 *
 * \param[in]	dim		dimension
 * \returns		number of objects of the dimension contained in the ref elem
 */
	static size_t substruct_num(int dim)
	{
	//	number of interior substructure Geometric Objects
		size_t vSubStructNum[MAXDIM+1];

		vSubStructNum[VERTEX] = 0;	// no additional vertices in the substructure
		vSubStructNum[EDGE] = 2;
		vSubStructNum[FACE] = 8;
		vSubStructNum[VOLUME] = 4;

		return vSubStructNum[dim];
	}

/// returns the number of objects of dim for a sub-geometric object of the implicit interior substructure
/**
 * The octahedral reference element contains an implicit interior
 * substructure that is constructed by several geometric objects, that
 * are mapped by a reference element by themselves. This method returns how
 * many (sub-)geometric objects of a given dimension are contained in the
 * (sub-)geometric object of the implicit interior substructure of
 * this reference element.
 *
 * \param[in]	dim_i		dimension of sub geometric object
 * \param[in]	i			number of sub geometric object
 * \param[in]	dim_j		dimension for elems contained in the sub-object
 * \returns		number of objects of the dimension dim_j that are
 * 				contained in the i*th (sub-)geom object of dimension dim_i
 */
	static size_t substruct_num(int dim_i, size_t i, int dim_j)
	{
	// 	number of interior substructure Geometric Objects
		size_t vSubStructSubNum[MAXDIM+1][MAXSUBSTRUCTOBJECTS][MAXDIM+1];

		for(size_t i = 0; i < substruct_num(VOLUME); ++i)
		{
			vSubStructSubNum[VOLUME][i][VERTEX] = 4;
			vSubStructSubNum[VOLUME][i][EDGE] = 6;
			vSubStructSubNum[VOLUME][i][FACE] = 4;
			vSubStructSubNum[VOLUME][i][VOLUME] = 1;
		}

		for(size_t i = 0; i < substruct_num(FACE); ++i)
		{
			vSubStructSubNum[FACE][i][VERTEX] = 3;
			vSubStructSubNum[FACE][i][EDGE] = 3;
			vSubStructSubNum[FACE][i][FACE] = 1;
			vSubStructSubNum[FACE][i][VOLUME] = 1;
		}

		for(size_t i = 0; i < substruct_num(EDGE); ++i)
		{
			vSubStructSubNum[EDGE][i][VERTEX] = 2;
			vSubStructSubNum[EDGE][i][EDGE] = 1;
			vSubStructSubNum[EDGE][i][FACE] = 2;
			vSubStructSubNum[EDGE][i][VOLUME] = 1;
		}

		for(size_t i = 0; i < substruct_num(VERTEX); ++i)
		{
			vSubStructSubNum[VERTEX][i][VERTEX] = 1;
			vSubStructSubNum[VERTEX][i][EDGE] = 3;
			vSubStructSubNum[VERTEX][i][FACE] = 3;
			vSubStructSubNum[VERTEX][i][VOLUME] = 1;
		}

		return vSubStructSubNum[dim_i][i][dim_j];
	}

/// returns the 'from' and 'to' corner indices for a scvf (overridden function from fv1_traits_most_common)
	static size_t scvf_from_to
	(
		const ReferenceOctahedron& refElem, ///< the reference element object
		size_t i, ///< index of the scvf
		size_t ft ///< 0 = from, 1 = to
	)
	{ // map according to the order defined in ComputeSCVFMidID
		static size_t from_to_ind [24][2] =
		{
			{1, 2},	// face 0
			{2, 1},	// face 1
			
			{2, 3},	// face 2
			{3, 2},	// face 3
			
			{3, 1},	// face 4
			{1, 3},	// face 5
			
			{1, 5},	// face 6
			{1, 0},	// face 7
			
			{2, 5},	// face 8
			{2, 0},	// face 9
			
			{3, 5},	// face 10
			{3, 0},	// face 11
			
			{1, 3},	// face 12
			{3, 1},	// face 13
			
			{3, 4},	// face 14
			{4, 3},	// face 15
			
			{4, 1},	// face 16
			{1, 4},	// face 17
			
			{1, 5},	// face 18
			{1, 0},	// face 19
			
			{3, 5},	// face 20
			{3, 0},	// face 21
			
			{4, 5},	// face 22
			{4, 0}	// face 23
		};
		return from_to_ind [i] [ft];
	}
	
///	returns the node id for a scv (overridden function from fv1_traits_most_common)
	static size_t scv_node_id
	(
		const ReferenceOctahedron& refElem, ///< the reference element object
		size_t i ///< index of the scv
	)
	{ // map according to order defined in ComputeSCVMidID
		static size_t node_id [16] =
		{
			1,	// scv 0
			1,	// scv 1
			
			2,	// scv 2
			2,	// scv 3
			
			3,	// scv 4
			3,	// scv 5
			
			5,	// scv 6
			0,	// scv 7
			
			1,	// scv 8
			1,	// scv 9
			
			3,	// scv 10
			3,	// scv 11
			
			4,	// scv 12
			4,	// scv 13
			
			5,	// scv 14
			0	// scv 15
		};
		return node_id [i];
	}
};
/// Octahedra: the FV1 traits
template <> struct fv1_traits<ReferenceOctahedron, 3> : public fv1_traits_ReferenceOctahedron {};

////////////////////////////////////////////////////////////////////////////////
// Dimension dependent traits DIM FV1
////////////////////////////////////////////////////////////////////////////////

///	Base for the Traits for Finite Volumes for a generic element of the fixed dimensionalities
template <int TDim, int TWorldDim> struct fv1_dim_traits_base
{
///	dimension of reference element
	static const int dim = TDim;

///	generic reference element type
	typedef DimReferenceElement<dim> ref_elem_type;
	
///	returns the number of the SCV
	static void dim_get_num_SCV_and_SCVF
	(
		const ref_elem_type& refElem, ///< reference element object
		ReferenceObjectID roid, ///< reference element object id
		size_t& numSCV, ///< to write the number of the SCVs
		size_t& numSCVF ///< to write the number of the SCVFs
	)
	{
		if(roid != ROID_PYRAMID && roid != ROID_OCTAHEDRON)
		{
			numSCV  = refElem.num(0);
			numSCVF = refElem.num(1);
		}
		else if(dim == 3 && roid == ROID_PYRAMID)
		{
			UG_WARNING("Pyramid Finite Volume Geometry for 1st order currently  "
					"implemented in DimFV1Geom EXPERIMENTATLLY. Please contact "
					"Martin Stepniewski or Andreas Vogel if you see this message.")

			numSCV  = 4*refElem.num(1);
			numSCVF = 2*refElem.num(1);
		}
		else if(dim == 3 && roid == ROID_OCTAHEDRON)
		{
		// 	Case octahedron
			numSCV  = 16;
			numSCVF = 24;
		}
		else
			UG_THROW ("fv1_dim_traits_base: Unsupported combination of dimension and reference element.");
	}
	
/// returns the 'from' and 'to' corner indices for a scvf
	static void get_dim_scvf_from_to
	(
		const ref_elem_type& refElem, ///< reference element object
		ReferenceObjectID roid, ///< reference element object id
		size_t i, ///< index of the scvf
		size_t& From, ///< to write the from-index
		size_t& To ///< to write the to-index
	)
	{
		if (roid != ROID_PYRAMID && roid != ROID_OCTAHEDRON)
		{
			From = refElem.id(1, i, 0, 0);
			To = refElem.id(1, i, 0, 1);
		}
		// special case pyramid (scvf not mappable by edges)
		else if (dim == 3 && roid == ROID_PYRAMID)
		{
		// 	map according to order defined in ComputeSCVFMidID
			From = fv1_traits_ReferencePyramid::scvf_from_to ((const ReferencePyramid&) refElem, i, 0);
			To = fv1_traits_ReferencePyramid::scvf_from_to ((const ReferencePyramid&) refElem, i, 1);
		}
		//	special case octahedron (scvf not mappable by edges)
		else if(dim == 3 && roid == ROID_OCTAHEDRON)
		{
		// 	map according to order defined in ComputeSCVFMidID
			From = fv1_traits_ReferenceOctahedron::scvf_from_to ((const ReferenceOctahedron&) refElem, i, 0);
			To = fv1_traits_ReferenceOctahedron::scvf_from_to ((const ReferenceOctahedron&) refElem, i, 1);
		}
		else
			UG_THROW ("fv1_dim_traits_base: Unsupported combination of dimension and reference element.");
	}
	
///	returns the node id for a scv
	static size_t dim_scv_node_id
	(
		const ref_elem_type& refElem, ///< reference element object
		ReferenceObjectID roid, ///< reference element object id
		size_t i ///< index of the scv
	)
	{
		// "classical" elements
		if (roid != ROID_PYRAMID && roid != ROID_OCTAHEDRON)
			return i;
		
		// special case pyramid (scv not mappable by corners)
		if(dim == 3 && roid == ROID_PYRAMID)
			return fv1_traits_ReferencePyramid::scv_node_id ((const ReferencePyramid&) refElem, i);
		
		// special case octahedron (scvf not mappable by edges)
		if(dim == 3 && roid == ROID_OCTAHEDRON)
			return fv1_traits_ReferenceOctahedron::scv_node_id ((const ReferenceOctahedron&) refElem, i);
		
		UG_THROW ("fv1_dim_traits_base: Unsupported combination of dimension and reference element.");
	}
};

///	Traits for Finite Volumes for a generic element of the fixed dimensionalities
template <int TDim, int TWorldDim> struct fv1_dim_traits;

template <> struct fv1_dim_traits<1, 1> : public fv1_traits<ReferenceEdge, 1>, public fv1_dim_traits_base<1, 1> {};
template <> struct fv1_dim_traits<1, 2> : public fv1_traits<ReferenceEdge, 2>, public fv1_dim_traits_base<1, 2> {};
template <> struct fv1_dim_traits<1, 3> : public fv1_traits<ReferenceEdge, 3>, public fv1_dim_traits_base<1, 3> {};

template <> struct fv1_dim_traits<2, 2> : public fv1_traits_ReferenceFace2d, public fv1_dim_traits_base<2, 2> {};
template <> struct fv1_dim_traits<2, 3> : public fv1_traits_ReferenceFace3d, public fv1_dim_traits_base<2, 3>
{
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

template <> struct fv1_dim_traits<3, 3>	: public fv1_traits_ReferenceVolume, public fv1_dim_traits_base<3, 3> {};

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

struct hfv1_traits_ReferenceVolume
{
	static void NormalOnSCVF(MathVector<3>& outNormal, const MathVector<3>* vCornerCoords)
	{ElementNormal<ReferenceTriangle, 3>(outNormal, vCornerCoords);}

	typedef ReferenceTetrahedron scv_type;
};

template <> struct hfv1_traits<ReferenceTetrahedron, 3> : public hfv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;
};

template <> struct hfv1_traits<ReferencePrism, 3> : public hfv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;
};

template <> struct hfv1_traits<ReferencePyramid, 3> : public hfv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 10;
};

template <> struct hfv1_traits<ReferenceHexahedron, 3> : public hfv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;
};

template <> struct hfv1_traits<ReferenceOctahedron, 3> : public hfv1_traits_ReferenceVolume
{
	const static size_t NumCornersOfSCVF = 3;
	const static size_t MaxNumCornersOfSCV = 8;
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

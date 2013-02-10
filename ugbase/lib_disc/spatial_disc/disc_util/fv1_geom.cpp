/*
 * fv1_geom.cpp
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */


#include "common/util/provider.h"
#include "fv1_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/quadrature/quadrature.h"
#include "lib_algebra/common/operations_vec.h"

namespace ug{


/**
 * \tparam	dim			dimension of coordinates
 * \tparam	TRefElem	Reference element type
 * \tparam	maxMid		Maximum number of elements for all dimensions
 */
template <int dim, typename TRefElem, int maxMid>
static void ComputeMidPoints(const TRefElem& rRefElem,
                             const MathVector<dim> vCorner[],
                             MathVector<dim> vvMid[][maxMid])
{
// 	compute local midpoints for all geometric objects with  0 < d <= dim
	for(int d = 1; d <= dim; ++d)
	{
	// 	loop geometric objects of dimension d
		for(size_t i = 0; i < rRefElem.num(d); ++i)
		{
		// 	set first node
			const size_t coID0 = rRefElem.id(d, i, 0, 0);
			vvMid[d][i] = vCorner[coID0];

		// 	add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = rRefElem.id(d, i, 0, j);
				vvMid[d][i] += vCorner[coID];
			}

		// 	scale for correct averaging
			vvMid[d][i] *= 1./(rRefElem.num(d, i, 0));
		}
	}

	// for PYRAMIDS: add midpoints of imaginary faces, edges and volumes
	// resulting from the division into two tetrahedra alongside x==y
	if (rRefElem.reference_object_id() == ROID_PYRAMID)
	{
		// diagonal 2->0, diagonal 0->2
		VecScaleAdd(vvMid[1][rRefElem.num(1)], 0.5, vCorner[2], 0.5, vCorner[0]);
		VecScaleAdd(vvMid[1][rRefElem.num(1)+1], 0.5, vCorner[0], 0.5, vCorner[2]);

		// subface 0,1,2; subface 0,2,3; face 0,4,2; face 0,2,4
		vvMid[2][rRefElem.num(2)] = vCorner[0];
		vvMid[2][rRefElem.num(2)] += vCorner[1];
		vvMid[2][rRefElem.num(2)] += vCorner[2];
		vvMid[2][rRefElem.num(2)] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+1] = vCorner[0];
		vvMid[2][rRefElem.num(2)+1] += vCorner[2];
		vvMid[2][rRefElem.num(2)+1] += vCorner[3];
		vvMid[2][rRefElem.num(2)+1] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+2] = vCorner[0];
		vvMid[2][rRefElem.num(2)+2] += vCorner[4];
		vvMid[2][rRefElem.num(2)+2] += vCorner[2];
		vvMid[2][rRefElem.num(2)+2] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+3] = vCorner[0];
		vvMid[2][rRefElem.num(2)+3] += vCorner[2];
		vvMid[2][rRefElem.num(2)+3] += vCorner[4];
		vvMid[2][rRefElem.num(2)+3] *= 1.0/3.0;

		// subvolume 0,1,2,4; subvolume 0,2,3,4

		vvMid[3][rRefElem.num(3)] = vCorner[0];
		vvMid[3][rRefElem.num(3)] += vCorner[1];
		vvMid[3][rRefElem.num(3)] += vCorner[2];
		vvMid[3][rRefElem.num(3)] += vCorner[4];
		vvMid[3][rRefElem.num(3)] *= 0.25;

		vvMid[3][rRefElem.num(3)+1] = vCorner[0];
		vvMid[3][rRefElem.num(3)+1] += vCorner[2];
		vvMid[3][rRefElem.num(3)+1] += vCorner[3];
		vvMid[3][rRefElem.num(3)+1] += vCorner[4];
		vvMid[3][rRefElem.num(3)+1] *= 0.25;
	}
}

/**
 * \param[in]	i		indicates that scvf corresponds to i'th edge of ref elem
 */
template <typename TRefElem>
static void ComputeSCVFMidID(const TRefElem& rRefElem,
                                   MidID vMidID[], int i)
{
	static const int dim = TRefElem::dim;

	if (rRefElem.reference_object_id() != ROID_PYRAMID)
	{
		//	set mid ids
		{
			// 	start at edge midpoint
			vMidID[0] = MidID(1,i);

			// 	loop up dimension
			if(dim == 2)
			{
				vMidID[1] = MidID(dim, 0); // center of element
			}
			else if (dim == 3)
			{
				vMidID[1] = MidID(2, rRefElem.id(1, i, 2, 0)); // side 0
				vMidID[2] = MidID(dim, 0); // center of element
				vMidID[3] = MidID(2, rRefElem.id(1, i, 2, 1)); // side 1
			}
		}
	}
	// pyramid here
	else
	{
		switch (i)
		{
			// scvf of edge 0
			case 0:	vMidID[0] = MidID(1,0);	// edge 0
					vMidID[1] = MidID(2,5);	// subface 0/0
					vMidID[2] = MidID(3,1);	// subvolume 0/0
					vMidID[3] = MidID(2,1); // face 1
					break;
			// scvf of edge 1
			case 1:	vMidID[0] = MidID(1,1);	// edge 1
					vMidID[1] = MidID(2,5);	// subface 0/0
					vMidID[2] = MidID(3,1);	// subvolume 0/0
					vMidID[3] = MidID(2,2); // face 2
					break;
			// scvf of diagonal 2->0
			case 2:	vMidID[0] = MidID(1,8);	// diagonal 2->0
					vMidID[1] = MidID(2,5);	// subface 0/0
					vMidID[2] = MidID(3,1);	// subvolume 0/0
					vMidID[3] = MidID(2,7); // face 0,4,2
					break;
			// scvf of edge 4 in subvolume 0/0
			case 3:	vMidID[0] = MidID(1,4);	// edge 4
					vMidID[1] = MidID(2,1); // face 1
					vMidID[2] = MidID(3,1);	// subvolume 0/0
					vMidID[3] = MidID(2,7); // face 0,4,2
					break;
			// scvf of edge 5
			case 4:	vMidID[0] = MidID(1,5);	// edge 5
					vMidID[1] = MidID(2,2);	// face 2
					vMidID[2] = MidID(3,1);	// subvolume 0/0
					vMidID[3] = MidID(2,1); // face 1
					break;
			// scvf of edge 6 in subvolume 0/0
			case 5:	vMidID[0] = MidID(1,6);	// edge 6
					vMidID[1] = MidID(2,7); // face 0,4,2
					vMidID[2] = MidID(3,1);	// subvolume 0/0
					vMidID[3] = MidID(2,2);	// face 2
					break;
			// edge 0->2
			case 6:	vMidID[0] = MidID(1,9);	// edge 0->2
					vMidID[1] = MidID(2,6);	// subface 1/0
					vMidID[2] = MidID(3,2);	// subvolume 1/0
					vMidID[3] = MidID(2,8);	// face 0,2,4
					break;
			// scvf of edge 2
			case 7:	vMidID[0] = MidID(1,2);	// edge 2
					vMidID[1] = MidID(2,6);	// subface 1/0
					vMidID[2] = MidID(3,2);	// subvolume 1/0
					vMidID[3] = MidID(2,3); // face 3
					break;
			// scvf of edge 3
			case 8:	vMidID[0] = MidID(1,3);	// edge 3
					vMidID[1] = MidID(2,6);	// subface 1/0
					vMidID[2] = MidID(3,2);	// subvolume 1/0
					vMidID[3] = MidID(2,4); // face 4
					break;
			// scvf of edge 4 in subvolume 1/0
			case 9:	vMidID[0] = MidID(1,4);	// edge 4
					vMidID[1] = MidID(2,8);	// face 0,2,4
					vMidID[2] = MidID(3,2);	// subvolume 1/0
					vMidID[3] = MidID(2,4); // face 4
					break;
			// scvf of edge 6 in subvolume 1/0
			case 10:vMidID[0] = MidID(1,6);	// edge 6
					vMidID[1] = MidID(2,3);	// face 3
					vMidID[2] = MidID(3,2);	// subvolume 1/0
					vMidID[1] = MidID(2,8);	// face 0,2,4
					break;
			// scvf of edge 7
			case 11:vMidID[0] = MidID(1,7);	// edge 7
					vMidID[3] = MidID(2,4); // face 4
					vMidID[2] = MidID(3,2);	// subvolume 1/0
					vMidID[1] = MidID(2,3);	// face 3
					break;
			default:UG_THROW("Pyramid only has 12 SCVFs (no. 0-11), but requested no. " << i << ".");
					break;
		}
	}
}

/**
 * \param[in]	i		indicates that scvf corresponds to i'th corner of ref elem
 */
template <typename TRefElem>
static void ComputeSCVMidID(const TRefElem& rRefElem,
                            MidID vMidID[], int i)
{
	static const int dim = TRefElem::dim;

	if (rRefElem.reference_object_id() != ROID_PYRAMID)
	{
		if(dim == 1)
		{
			vMidID[0] = MidID(0, i); // set node as corner of scv
			vMidID[1] = MidID(dim, 0);	// center of element
		}
		else if(dim == 2)
		{
			vMidID[0] = MidID(0, i); // set node as corner of scv
			vMidID[1] = MidID(1, rRefElem.id(0, i, 1, 0)); // edge 1
			vMidID[2] = MidID(dim, 0);	// center of element
			vMidID[3] = MidID(1, rRefElem.id(0, i, 1, 1)); // edge 2
		}
		else if(dim == 3)
		{
			vMidID[0] = MidID(0, i); // set node as corner of scv
			vMidID[1] = MidID(1, rRefElem.id(0, i, 1, 1)); // edge 1
			vMidID[2] = MidID(2, rRefElem.id(0, i, 2, 0)); // face 0
			vMidID[3] = MidID(1, rRefElem.id(0, i, 1, 0)); // edge 0
			vMidID[4] = MidID(1, rRefElem.id(0, i, 1, 2)); // edge 2
			vMidID[5] = MidID(2, rRefElem.id(0, i, 2, 2)); // face 2
			vMidID[6] = MidID(dim, 0);	// center of element
			vMidID[7] = MidID(2, rRefElem.id(0, i, 2, 1)); // face 1
		}
		else {UG_THROW("Dimension higher that 3 not implemented.");}
	}
	// pyramid here
	else
	{
		switch (i)
		{
			// scv of corner 0 in subvolume 0/0
			case 0:	vMidID[0] = MidID(0,0);	// corner 0
					vMidID[1] = MidID(1,0);	// edge 0
					vMidID[2] = MidID(2,5);	// subface 0/0
					vMidID[3] = MidID(1,8); // edge 2->0
					vMidID[4] = MidID(1,4);	// edge 4
					vMidID[5] = MidID(2,1); // face 1
					vMidID[6] = MidID(3,1);	// subvolume 0/0
					vMidID[7] = MidID(2,7); // face 0,4,2
					break;
			// scv of corner 1
			case 1:	vMidID[0] = MidID(0,1);	// corner 1
					vMidID[1] = MidID(1,1);	// edge 1
					vMidID[2] = MidID(2,5);	// subface 0/0
					vMidID[3] = MidID(1,0);	// edge 0
					vMidID[4] = MidID(1,5);	// edge 5
					vMidID[5] = MidID(2,2);	// face 2
					vMidID[6] = MidID(3,1);	// subvolume 0/0
					vMidID[7] = MidID(2,1); // face 1
					break;
			// scv of corner 2 in subvolume 0/0
			case 2:	vMidID[0] = MidID(0,2);	// corner 2
					vMidID[1] = MidID(1,8); // edge 2->0
					vMidID[2] = MidID(2,5);	// subface 0/0
					vMidID[3] = MidID(1,1);	// edge 1
					vMidID[4] = MidID(1,6);	// edge 6
					vMidID[5] = MidID(2,7); // face 0,4,2
					vMidID[6] = MidID(3,1);	// subvolume 0/0
					vMidID[7] = MidID(2,2);	// face 2
					break;
			// scv of corner 4 in subvolume 0/0
			case 3:	vMidID[0] = MidID(0,4);	// corner 4
					vMidID[1] = MidID(1,5);	// edge 5
					vMidID[2] = MidID(2,1); // face 1
					vMidID[3] = MidID(1,4); // edge 4
					vMidID[4] = MidID(1,6);	// edge 6
					vMidID[5] = MidID(2,2); // face 2
					vMidID[6] = MidID(3,1);	// subvolume 0/0
					vMidID[7] = MidID(2,7); // face 0,4,2
					break;
			// scv of corner 0 in subvolume 1/0
			case 4:	vMidID[0] = MidID(0,0);	// corner 0
					vMidID[1] = MidID(1,9);	// edge 0->2
					vMidID[2] = MidID(2,6);	// subface 1/0
					vMidID[3] = MidID(1,3); // edge 3
					vMidID[4] = MidID(1,4);	// edge 4
					vMidID[5] = MidID(2,8); // face 0,2,4
					vMidID[6] = MidID(3,2);	// subvolume 1/0
					vMidID[7] = MidID(2,4); // face 4
					break;
			// scv of corner 2 in subvolume 1/0
			case 5:	vMidID[0] = MidID(0,2);	// corner 2
					vMidID[1] = MidID(1,2);	// edge 2
					vMidID[2] = MidID(2,6);	// subface 1/0
					vMidID[3] = MidID(1,9); // edge 0->2
					vMidID[4] = MidID(1,6);	// edge 6
					vMidID[5] = MidID(2,3); // face 3
					vMidID[6] = MidID(3,2);	// subvolume 1/0
					vMidID[7] = MidID(2,8); // face 0,2,4
					break;
			// scv of corner 3
			case 6:	vMidID[0] = MidID(0,3);	// corner 3
					vMidID[1] = MidID(1,3);	// edge 3
					vMidID[2] = MidID(2,6);	// subface 1/0
					vMidID[3] = MidID(1,2); // edge 2
					vMidID[4] = MidID(1,7);	// edge 7
					vMidID[5] = MidID(2,4); // face 4
					vMidID[6] = MidID(3,2);	// subvolume 1/0
					vMidID[7] = MidID(2,3); // face 3
					break;
			// scv of corner 4 in subvolume 1/0
			case 7:	vMidID[0] = MidID(0,4);	// corner 4
					vMidID[1] = MidID(1,6);	// edge 6
					vMidID[2] = MidID(2,8); // face 0,2,4
					vMidID[3] = MidID(1,4); // edge 4
					vMidID[4] = MidID(1,7);	// edge 7
					vMidID[5] = MidID(2,3); // face 3
					vMidID[6] = MidID(3,2);	// subvolume 1/0
					vMidID[7] = MidID(2,4); // face 4
					break;
			default:UG_THROW("Pyramid only has 8 SCVs (no. 0-7), but requested no. " << i << ".");
					break;
		}
	}
}

/**
 * \param[in]	i		indicates that scvf corresponds to i'th corner of ref elem
 */
template <typename TRefElem>
static void ComputeBFMidID(const TRefElem& rRefElem, int side,
                            MidID vMidID[], int co)
{
	static const int dim = TRefElem::dim;

	if (rRefElem.reference_object_id() != ROID_PYRAMID || side != 0)
	{
		//	number of corners of side
		const int coOfSide = rRefElem.num(dim-1, side, 0);

		// 	set mid ids
		if(dim == 2)
		{
			vMidID[co%2] = MidID(0, rRefElem.id(1, side, 0, co)); // corner of side
			vMidID[(co+1)%2] = MidID(1, side); // side midpoint
		}
		else if (dim == 3)
		{
			vMidID[0] = MidID(0, rRefElem.id(2, side, 0, co)); // corner of side
			vMidID[1] = MidID(1, rRefElem.id(2, side, 1, co)); // edge co
			vMidID[2] = MidID(2, side); // side midpoint
			vMidID[3] = MidID(1, rRefElem.id(2, side, 1, (co -1 + coOfSide)%coOfSide)); // edge co-1
		}
	}
	// bottom side of pyramid here
	else
	{
		switch (co)
		{
			// bf of corner 0 in subface 0/0
			case 0:	vMidID[0] = MidID(0,0);	// corner 0
					vMidID[1] = MidID(1,8); // edge 2->0
					vMidID[2] = MidID(2,5);	// subface 0/0
					vMidID[3] = MidID(1,0);	// edge 0
					break;
			// bf of corner 1
			case 1:	vMidID[0] = MidID(0,1);	// corner 1
					vMidID[1] = MidID(1,0);	// edge 0
					vMidID[2] = MidID(2,5);	// subface 0/0
					vMidID[3] = MidID(1,1);	// edge 1
					break;
			// bf of corner 2 in subvolume 0/0
			case 2:	vMidID[0] = MidID(0,2);	// corner 2
					vMidID[1] = MidID(1,1);	// edge 1
					vMidID[2] = MidID(2,5);	// subface 0/0
					vMidID[3] = MidID(1,8); // edge 2->0
					break;
			// bf of corner 0 in subvolume 1/0
			case 3:	vMidID[0] = MidID(0,0);	// corner 0
					vMidID[1] = MidID(1,3); // edge 3
					vMidID[2] = MidID(2,6);	// subface 1/0
					vMidID[3] = MidID(1,9);	// edge 0->2
					break;
			// bf of corner 2 in subvolume 1/0
			case 4:	vMidID[0] = MidID(0,2);	// corner 2
					vMidID[1] = MidID(1,9); // edge 0->2
					vMidID[2] = MidID(2,6);	// subface 1/0
					vMidID[3] = MidID(1,2);	// edge 2
					break;
			// bf of corner 3
			case 5:	vMidID[0] = MidID(0,3);	// corner 3
					vMidID[1] = MidID(1,2);	// edge 2
					vMidID[2] = MidID(2,6);	// subface 1/0
					vMidID[3] = MidID(1,3); // edge 3
					break;
			default:UG_THROW("Pyramid only has 6 BFs on bottom side (no. 0-5), but requested no. " << co << ".");
					break;
		}
	}
}

template <int dim, int maxMid>
static void CopyCornerByMidID(MathVector<dim> vCorner[],
                              const MidID vMidID[],
                              MathVector<dim> vvMidPos[][maxMid],
                              const size_t numCo)
{
	for(size_t i = 0; i < numCo; ++i)
	{
		const size_t d = vMidID[i].dim;
		const size_t id = vMidID[i].id;
		vCorner[i] = vvMidPos[d][id];
	}
}

////////////////////////////////////////////////////////////////////////////////
// FV1 Geometry for Reference Element Type
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
FV1Geometry<TElem, TWorldDim>::
FV1Geometry()
	: m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()),
	  m_rTrialSpace(Provider<local_shape_fct_set_type>::get())
{
	update_local_data();
}

template <typename TElem, int TWorldDim>
void FV1Geometry<TElem, TWorldDim>::
update_local_data()
{
// 	set corners of element as local centers of nodes
	for(size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_vvLocMid[0][i] = m_rRefElem.corner(i);

//	compute local midpoints
	ComputeMidPoints<dim, ref_elem_type, maxMid>(m_rRefElem, m_vvLocMid[0], m_vvLocMid);

// 	set up local information for SubControlVolumeFaces (scvf)
	for(size_t i = 0; i < num_scvf(); ++i)
	{

	//	this scvf separates the given nodes
		if (m_rRefElem.REFERENCE_OBJECT_ID != ROID_PYRAMID)
		{
			m_vSCVF[i].From = m_rRefElem.id(1, i, 0, 0);
			m_vSCVF[i].To = m_rRefElem.id(1, i, 0, 1);
		}
		// special case pyramid (scvf not mappable by edges)
		else
		{
			// map according to order defined in ComputeSCVFMidID
			m_vSCVF[i].From = ((i>6 && i%3) ? (i%3)+1 : i%3);
			m_vSCVF[i].To = i%6 > 2 ? 4 : ((i+1)%3 + (i>5 && i<8 ? 1 : 0));
		}

	//	compute mid ids of the scvf
		ComputeSCVFMidID(m_rRefElem, m_vSCVF[i].vMidID, i);

	// 	copy local corners of scvf
		CopyCornerByMidID<dim, maxMid>(m_vSCVF[i].vLocPos, m_vSCVF[i].vMidID, m_vvLocMid, SCVF::numCo);

	// 	integration point
		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].vLocPos, SCVF::numCo);
	}

// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the element
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	store associated node
		if (m_rRefElem.REFERENCE_OBJECT_ID != ROID_PYRAMID)
		{
			m_vSCV[i].nodeId = i;
		}
		// special case pyramid (scv not mappable by corners)
		else
		{
			// map according to order defined in ComputeSCVMidID
			m_vSCV[i].nodeId = i<3 ? i : (i<5 ? (i+1)%5 : i-3);
		}

	//	compute mid ids scv
		ComputeSCVMidID(m_rRefElem, m_vSCV[i].midId, i);

	// 	copy local corners of scv
		CopyCornerByMidID<dim, maxMid>(m_vSCV[i].vLocPos, m_vSCV[i].midId, m_vvLocMid, m_vSCV[i].num_corners());
	}


// 	compute Shapes and Derivatives
	for(size_t i = 0; i < num_scvf(); ++i)
	{
		m_rTrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].local_ip());
		m_rTrialSpace.grads(&(m_vSCVF[i].vLocalGrad[0]), m_vSCVF[i].local_ip());
	}

	for(size_t i = 0; i < num_scv(); ++i)
	{
		m_rTrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].local_ip());
		m_rTrialSpace.grads(&(m_vSCV[i].vLocalGrad[0]), m_vSCV[i].local_ip());
	}

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();
}

/// update data for given element
template <typename TElem, int TWorldDim>
void FV1Geometry<TElem, TWorldDim>::
update(TElem* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
// 	if already update for this element, do nothing
	if(m_pElem == elem) return; else m_pElem = elem;

// 	remember global position of nodes
	for(size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_vvGloMid[0][i] = vCornerCoords[i];

//	compute global midpoints
	ComputeMidPoints<worldDim, ref_elem_type, maxMid>(m_rRefElem, m_vvGloMid[0], m_vvGloMid);

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	// 	copy local corners of scvf
		CopyCornerByMidID<worldDim, maxMid>(m_vSCVF[i].vGloPos, m_vSCVF[i].vMidID, m_vvGloMid, SCVF::numCo);

	// 	integration point
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].vGloPos, SCVF::numCo);

	// 	normal on scvf
		traits::NormalOnSCVF(m_vSCVF[i].Normal, m_vSCVF[i].vGloPos, m_vvGloMid[0]);
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
	// 	copy global corners
		CopyCornerByMidID<worldDim, maxMid>(m_vSCV[i].vGloPos, m_vSCV[i].midId, m_vvGloMid, m_vSCV[i].num_corners());

	// 	compute volume of scv
		m_vSCV[i].Vol = ElementSize<scv_type, worldDim>(m_vSCV[i].vGloPos);
	}

// 	Shapes and Derivatives
	m_mapping.update(vCornerCoords);

//	if mapping is linear, compute jacobian only once and copy
	if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
	{
		MathMatrix<worldDim,dim> JtInv;
		m_mapping.jacobian_transposed_inverse(JtInv, m_vSCVF[0].local_ip());
		const number detJ = m_mapping.sqrt_gram_det(m_vSCVF[0].local_ip());

		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_vSCVF[i].JtInv = JtInv;
			m_vSCVF[i].detj = detJ;
		}

		for(size_t i = 0; i < num_scv(); ++i)
		{
			m_vSCV[i].JtInv = JtInv;
			m_vSCV[i].detj = detJ;
		}
	}
//	else compute jacobian for each integration point
	else
	{
		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_mapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].local_ip());
			m_vSCVF[i].detj = m_mapping.sqrt_gram_det(m_vSCVF[i].local_ip());
		}
		for(size_t i = 0; i < num_scv(); ++i)
		{
			m_mapping.jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCV[i].local_ip());
			m_vSCV[i].detj = m_mapping.sqrt_gram_det(m_vSCV[i].local_ip());
		}
	}

//	compute global gradients
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t sh = 0 ; sh < scvf(i).num_sh(); ++sh)
			MatVecMult(m_vSCVF[i].vGlobalGrad[sh], m_vSCVF[i].JtInv, m_vSCVF[i].vLocalGrad[sh]);

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t sh = 0 ; sh < scv(i).num_sh(); ++sh)
			MatVecMult(m_vSCV[i].vGlobalGrad[sh], m_vSCV[i].JtInv, m_vSCV[i].vLocalGrad[sh]);

// 	Copy ip pos in list for SCVF
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vGlobSCVF_IP[i] = scvf(i).global_ip();

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(elem, vCornerCoords, ish);
}

template <typename TElem, int TWorldDim>
void FV1Geometry<TElem, TWorldDim>::
update_boundary_faces(TElem* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<VertexBase*> vVertex;
		CollectVertices(vVertex, grid, elem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<EdgeBase*> vEdges;
		CollectEdgesSorted(vEdges, grid, elem);
		vSubsetIndex.resize(vEdges.size());
		for(size_t i = 0; i < vEdges.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vEdges[i]);
	}
	if(dim == 3) {
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, elem);
		vSubsetIndex.resize(vFaces.size());
		for(size_t i = 0; i < vFaces.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vFaces[i]);
	}

//	loop requested subset
	typename std::map<int, std::vector<BF> >::iterator it;
	for (it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); ++it)
	{
	//	get subset index
		const int bndIndex = (*it).first;

	//	get vector of BF for element
		std::vector<BF>& vBF = (*it).second;

	//	clear vector
		vBF.clear();

	//	current number of bf
		size_t curr_bf = 0;

	//	loop sides of element
		for(size_t side = 0; side < vSubsetIndex.size(); ++side)
		{
		//	skip non boundary sides
			if(vSubsetIndex[side] != bndIndex) continue;

		//	number of corners of side (special case bottom side pyramid)
			const int coOfSide = (m_rRefElem.REFERENCE_OBJECT_ID != ROID_PYRAMID || side != 0)
								? m_rRefElem.num(dim-1, side, 0) : m_rRefElem.num(dim-1, side, 0) + 2;

		//	resize vector
			vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				BF& bf = vBF[curr_bf];

			//	set node id == scv this bf belongs to
				if (m_rRefElem.REFERENCE_OBJECT_ID != ROID_PYRAMID || side != 0)
					bf.nodeId = m_rRefElem.id(dim-1, side, 0, co);
				else
				{
					// map according to order defined in ComputeBFMidID
					bf.nodeId = m_rRefElem.id(dim-1, side, 0, (co % 3) + (co>3 ? 1 : 0));
				}

			//	Compute MidID for BF
				ComputeBFMidID(m_rRefElem, side, bf.vMidID, co);

			// 	copy corners of bf
				CopyCornerByMidID<dim, maxMid>(bf.vLocPos, bf.vMidID, m_vvLocMid, BF::numCo);
				CopyCornerByMidID<worldDim, maxMid>(bf.vGloPos, bf.vMidID, m_vvGloMid, BF::numCo);

			// 	integration point
				AveragePositions(bf.localIP, bf.vLocPos, BF::numCo);
				AveragePositions(bf.globalIP, bf.vGloPos, BF::numCo);

			// 	normal on scvf
				traits::NormalOnSCVF(bf.Normal, bf.vGloPos, m_vvGloMid[0]);

			//	compute volume
				bf.Vol = VecTwoNorm(bf.Normal);

				m_rTrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
				m_rTrialSpace.grads(&(bf.vLocalGrad[0]), bf.localIP);

				m_mapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
				bf.detj = m_mapping.sqrt_gram_det(bf.localIP);

				for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
					MatVecMult(bf.vGlobalGrad[sh], bf.JtInv, bf.vLocalGrad[sh]);

			//	increase curr_bf
				++curr_bf;
			}
		}
	}
}

////////////////////////////////////////////////////////////////////////////////
// Dim-dependent Finite Volume Geometry
////////////////////////////////////////////////////////////////////////////////

template <int TDim, int TWorldDim>
void DimFV1Geometry<TDim, TWorldDim>::
update_local_data()
{
//	get reference element
	try{
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

// 	set corners of element as local centers of nodes
	for(size_t i = 0; i < rRefElem.num(0); ++i)
		m_vvLocMid[0][i] = rRefElem.corner(i);

//	compute local midpoints
	ComputeMidPoints<dim, DimReferenceElement<dim>, maxMid>
										(rRefElem, m_vvLocMid[0], m_vvLocMid);

//	set number of scvf / scv of this roid
	m_numSCV = (m_roid != ROID_PYRAMID) ? rRefElem.num(0) : 8; // number of corners
	m_numSCVF = (m_roid != ROID_PYRAMID) ? rRefElem.num(1) : 12; // number of edges

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one edge of the element
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	this scvf separates the given nodes
		if (m_roid != ROID_PYRAMID)
		{
			m_vSCVF[i].From = rRefElem.id(1, i, 0, 0);
			m_vSCVF[i].To = rRefElem.id(1, i, 0, 1);
		}
		// special case pyramid (scvf not mappable by edges)
		else
		{
			// map according to order defined in ComputeSCVFMidID
			m_vSCVF[i].From = ((i>6 && i%3) ? (i%3)+1 : i%3);
			m_vSCVF[i].To = i%6 > 2 ? 4 : ((i+1)%3 + (i>5 && i<8 ? 1 : 0));
		}


	//	compute mid ids of the scvf
		ComputeSCVFMidID(rRefElem, m_vSCVF[i].vMidID, i);

	// 	copy local corners of scvf
		CopyCornerByMidID<dim, maxMid>(m_vSCVF[i].vLocPos, m_vSCVF[i].vMidID, m_vvLocMid, SCVF::numCo);

	// 	integration point
		AveragePositions(m_vSCVF[i].localIP, m_vSCVF[i].vLocPos, SCVF::numCo);
	}

// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the element
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	store associated node
		if (m_roid != ROID_PYRAMID)
		{
			m_vSCV[i].nodeId = i;
		}
		// special case pyramid (scv not mappable by corners)
		else
		{
			// map according to order defined in ComputeSCVMidID
			m_vSCV[i].nodeId = i<3 ? i : (i<5 ? (i+1)%5 : i-3);
		}

	//	compute mid ids scv
		ComputeSCVMidID(rRefElem, m_vSCV[i].vMidID, i);

	// 	copy local corners of scv
		CopyCornerByMidID<dim, maxMid>(m_vSCV[i].vLocPos, m_vSCV[i].vMidID, m_vvLocMid, m_vSCV[i].num_corners());
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

	const LocalShapeFunctionSet<dim>& TrialSpace =
		LocalShapeFunctionSetProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, 1));

	m_nsh = TrialSpace.num_sh();

	for(size_t i = 0; i < num_scvf(); ++i)
	{
		m_vSCVF[i].numSH = TrialSpace.num_sh();
		TrialSpace.shapes(&(m_vSCVF[i].vShape[0]), m_vSCVF[i].localIP);
		TrialSpace.grads(&(m_vSCVF[i].vLocalGrad[0]), m_vSCVF[i].localIP);
	}

	for(size_t i = 0; i < num_scv(); ++i)
	{
		m_vSCV[i].numSH = TrialSpace.num_sh();
		TrialSpace.shapes(&(m_vSCV[i].vShape[0]), m_vSCV[i].vLocPos[0]);
		TrialSpace.grads(&(m_vSCV[i].vLocalGrad[0]), m_vSCV[i].vLocPos[0]);
	}

	}
	UG_CATCH_THROW("DimFV1Geometry: update failed.");

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();
}


/// update data for given element
template <int TDim, int TWorldDim>
void DimFV1Geometry<TDim, TWorldDim>::
update(GeometricObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

//	refresh local data, if different roid given
	if(m_roid != pElem->reference_object_id())
	{
	//	remember new roid
		m_roid = (ReferenceObjectID) pElem->reference_object_id();

	//	update local data
		update_local_data();
	}

//	get reference element
	try{
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

// 	remember global position of nodes
	for(size_t i = 0; i < rRefElem.num(0); ++i)
		m_vvGloMid[0][i] = vCornerCoords[i];

//	compute local midpoints
	ComputeMidPoints<worldDim, DimReferenceElement<dim>, maxMid>(rRefElem, m_vvGloMid[0], m_vvGloMid);

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	// 	copy local corners of scvf
		CopyCornerByMidID<worldDim, maxMid>(m_vSCVF[i].vGloPos, m_vSCVF[i].vMidID, m_vvGloMid, SCVF::numCo);

	// 	integration point
		AveragePositions(m_vSCVF[i].globalIP, m_vSCVF[i].vGloPos, SCVF::numCo);

	// 	normal on scvf
		traits::NormalOnSCVF(m_vSCVF[i].Normal, m_vSCVF[i].vGloPos, m_vvGloMid[0]);
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
	// 	copy global corners
		CopyCornerByMidID<worldDim, maxMid>(m_vSCV[i].vGloPos, m_vSCV[i].vMidID, m_vvGloMid, m_vSCV[i].num_corners());

	// 	compute volume of scv
		m_vSCV[i].Vol = ElementSize<scv_type, worldDim>(m_vSCV[i].vGloPos);
	}

//	get reference mapping
	DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
	rMapping.update(vCornerCoords);

	//\todo compute with on virt. call
//	compute jacobian for linear mapping
	if(rMapping.is_linear())
	{
		MathMatrix<worldDim,dim> JtInv;
		rMapping.jacobian_transposed_inverse(JtInv, m_vSCVF[0].local_ip());
		const number detJ = rMapping.sqrt_gram_det(m_vSCVF[0].local_ip());

		for(size_t i = 0; i < num_scvf(); ++i)
		{
			m_vSCVF[i].JtInv = JtInv;
			m_vSCVF[i].detj = detJ;
		}

		for(size_t i = 0; i < num_scv(); ++i)
		{
			m_vSCV[i].JtInv = JtInv;
			m_vSCV[i].detj = detJ;
		}
	}
//	else compute jacobian for each integration point
	else
	{
		for(size_t i = 0; i < num_scvf(); ++i)
		{
			rMapping.jacobian_transposed_inverse(m_vSCVF[i].JtInv, m_vSCVF[i].local_ip());
			m_vSCVF[i].detj = rMapping.sqrt_gram_det(m_vSCVF[i].local_ip());
		}
		for(size_t i = 0; i < num_scv(); ++i)
		{
			rMapping.jacobian_transposed_inverse(m_vSCV[i].JtInv, m_vSCV[i].local_ip());
			m_vSCV[i].detj = rMapping.sqrt_gram_det(m_vSCV[i].local_ip());
		}
	}

//	compute global gradients
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t sh = 0; sh < scvf(i).num_sh(); ++sh)
			MatVecMult(m_vSCVF[i].vGlobalGrad[sh], m_vSCVF[i].JtInv, m_vSCVF[i].vLocalGrad[sh]);

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t sh = 0; sh < scv(i).num_sh(); ++sh)
			MatVecMult(m_vSCV[i].vGlobalGrad[sh], m_vSCV[i].JtInv, m_vSCV[i].vLocalGrad[sh]);

// 	copy ip points in list (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vGlobSCVF_IP[i] = scvf(i).global_ip();

	}
	UG_CATCH_THROW("DimFV1Geometry: update failed.");

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}

template <int TDim, int TWorldDim>
void DimFV1Geometry<TDim, TWorldDim>::
update_boundary_faces(GeometricObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<VertexBase*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<EdgeBase*> vEdges;
		CollectEdgesSorted(vEdges, grid, pElem);
		vSubsetIndex.resize(vEdges.size());
		for(size_t i = 0; i < vEdges.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vEdges[i]);
	}
	if(dim == 3) {
		std::vector<Face*> vFaces;
		CollectFacesSorted(vFaces, grid, pElem);
		vSubsetIndex.resize(vFaces.size());
		for(size_t i = 0; i < vFaces.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vFaces[i]);
	}

	try{
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

	DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
	rMapping.update(vCornerCoords);

	const LocalShapeFunctionSet<dim>& TrialSpace =
		LocalShapeFunctionSetProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, 1));

//	loop requested subset
	typename std::map<int, std::vector<BF> >::iterator it;
	for (it=m_mapVectorBF.begin() ; it != m_mapVectorBF.end(); ++it)
	{
	//	get subset index
		const int bndIndex = (*it).first;

	//	get vector of BF for element
		std::vector<BF>& vBF = (*it).second;

	//	clear vector
		vBF.clear();

	//	current number of bf
		size_t curr_bf = 0;

	//	loop sides of element
		for(size_t side = 0; side < vSubsetIndex.size(); ++side)
		{
		//	skip non boundary sides
			if(vSubsetIndex[side] != bndIndex) continue;

		//	number of corners of side (special case bottom side pyramid)
			const int coOfSide = (pElem->reference_object_id() != ROID_PYRAMID || side != 0)
								? rRefElem.num(dim-1, side, 0) : rRefElem.num(dim-1, side, 0) + 2;
		//	resize vector
			vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				BF& bf = vBF[curr_bf];

			//	set node id == scv this bf belongs to
				if (pElem->reference_object_id() != ROID_PYRAMID || side != 0)
					bf.nodeId = rRefElem.id(dim-1, side, 0, co);
				else
				{
					// map according to order defined in ComputeBFMidID
					bf.nodeId = rRefElem.id(dim-1, side, 0, (co % 3) + (co>3 ? 1 : 0));
				}

			//	Compute MidID for BF
				ComputeBFMidID(rRefElem, side, bf.vMidID, co);

			// 	copy corners of bf
				CopyCornerByMidID<dim, maxMid>(bf.vLocPos, bf.vMidID, m_vvLocMid, BF::numCo);
				CopyCornerByMidID<worldDim, maxMid>(bf.vGloPos, bf.vMidID, m_vvGloMid, BF::numCo);

			// 	integration point
				AveragePositions(bf.localIP, bf.vLocPos, BF::numCo);
				AveragePositions(bf.globalIP, bf.vGloPos, BF::numCo);

			// 	normal on scvf
				traits::NormalOnSCVF(bf.Normal, bf.vGloPos, m_vvGloMid[0]);

			//	compute volume
				bf.Vol = VecTwoNorm(bf.Normal);

			//	compute shapes and grads
				bf.numSH = TrialSpace.num_sh();
				TrialSpace.shapes(&(bf.vShape[0]), bf.localIP);
				TrialSpace.grads(&(bf.vLocalGrad[0]), bf.localIP);

			//	get reference mapping
				rMapping.jacobian_transposed_inverse(bf.JtInv, bf.localIP);
				bf.detj = rMapping.sqrt_gram_det(bf.localIP);

			//	compute global gradients
				for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
					MatVecMult(bf.vGlobalGrad[sh], bf.JtInv, bf.vLocalGrad[sh]);

			//	increase curr_bf
				++curr_bf;
			}
		}
	}

	}
	UG_CATCH_THROW("DimFV1Geometry: update failed.");
}

////////////////////////////////////////////////////////////////////////////////
// FV1ManifoldBoundary
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
FV1ManifoldBoundary<TElem, TWorldDim>::
FV1ManifoldBoundary() : m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()), m_ssi(-1)
{
	// set corners of element as local centers of nodes
	for (size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_locMid[0][i] = m_rRefElem.corner(i);

	// compute local midpoints for all geometric objects with  0 < d <= dim
	for (int d = 1; d <= dim; ++d)
	{
		// loop geometric objects of dimension d
		for(size_t i = 0; i < m_rRefElem.num(d); ++i)
		{
			// set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_locMid[d][i] = m_locMid[0][coID0];

			// add corner coordinates of the corners of the geometric object
			for(size_t j = 1; j < m_rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_locMid[d][i] += m_locMid[0][coID];
			}

			// scale for correct averaging
			m_locMid[d][i] *= 1./(m_rRefElem.num(d, i, 0));
		}
	}

	// set up local information for Boundary Faces (bf)
	// each bf is associated to one corner of the element
	for (size_t i = 0; i < num_bf(); ++i)
	{
		m_vBF[i].nodeId = i;

		if (dim == 1) // Edge
		{
			m_vBF[i].midId[0] = MidID(0, i);	// set node as corner of bf
			m_vBF[i].midId[1] = MidID(dim, 0);	// center of bnd element
			
			// copy local corners of bf
			copy_local_corners(m_vBF[i]);
			
			// local integration point
			AveragePositions(m_vBF[i].localIP, m_vBF[i].vLocPos, 2);
		}
		else if (dim == 2)	// Quadrilateral
		{
			m_vBF[i].midId[0] = MidID(0, i); // set node as corner of bf
			m_vBF[i].midId[1] = MidID(1, m_rRefElem.id(0, i, 1, 0)); // edge 1
			m_vBF[i].midId[2] = MidID(dim, 0);	// center of bnd element
			m_vBF[i].midId[3] = MidID(1, m_rRefElem.id(0, i, 1, 1)); // edge 2
			
			// copy local corners of bf
			copy_local_corners(m_vBF[i]);
			
			// local integration point
			AveragePositions(m_vBF[i].localIP, m_vBF[i].vLocPos, 4);
		}
		else {UG_ASSERT(0, "Dimension higher that 2 not implemented.");}
	}

	/////////////
	// Shapes
	/////////////
	// A word of warning: This is only meaningful,
	// if the trial space is piecewise linear on tetrahedrons/triangles!
	for (size_t i = 0; i < num_bf(); ++i)
	{
		const LocalShapeFunctionSet<ref_elem_type::dim>& TrialSpace =
				LocalShapeFunctionSetProvider::
					get<ref_elem_type::dim>
					(ref_elem_type::REFERENCE_OBJECT_ID,
					 LFEID(LFEID::LAGRANGE, 1));

		const size_t num_sh = m_numBF;
		m_vBF[i].vShape.resize(num_sh);

		TrialSpace.shapes(&(m_vBF[i].vShape[0]), m_vBF[i].localIP);
	}

	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

	// 	loop Boundary Faces (BF)
	m_vLocBFIP.clear();
	for (size_t i = 0; i < num_bf(); ++i)
	{
	//	get current BF
		const BF& rBF = bf(i);

	// 	loop ips
		for (size_t ip = 0; ip < rBF.num_ip(); ++ip)
		{
			m_vLocBFIP.push_back(rBF.local_ip(ip));
		}
	}
}


/// update data for given element
template <typename TElem, int TWorldDim>
void FV1ManifoldBoundary<TElem, TWorldDim>::
update(TElem* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	// 	if already update for this element, do nothing
	if (m_pElem == elem) return;
	else m_pElem = elem;

	//	store subset index for geometry
	Grid& grid = *(ish->grid());
	if (worldDim == 2)
	{
		std::vector<EdgeBase*> vEdge;
		CollectEdgesSorted(vEdge, grid, elem);
		UG_ASSERT(vEdge.size(),"No edge contained in 1D manifold element!");
		m_ssi = ish->get_subset_index(vEdge[0]);
	}
	if (worldDim == 3)
	{
		std::vector<Face*> vFace;
		CollectFacesSorted(vFace, grid, elem);
		UG_ASSERT(vFace.size(),"No face contained in 2D manifold element!");
		m_ssi = ish->get_subset_index(vFace[0]);
	}

	// 	remember global position of nodes
	for (size_t i = 0; i < m_rRefElem.num(0); ++i)
		m_gloMid[0][i] = vCornerCoords[i];

	// 	compute global midpoints for all the other geometric objects (with  0 < d <= dim)
	for (int d = 1; d <= dim; ++d)
	{
		// 	loop geometric objects of dimension d
		for (size_t i = 0; i < m_rRefElem.num(d); ++i)
		{
			// set first node
			const size_t coID0 = m_rRefElem.id(d, i, 0, 0);
			m_gloMid[d][i] = m_gloMid[0][coID0];

		// 	add corner coordinates of the corners of the geometric object
			for (size_t j = 1; j < m_rRefElem.num(d, i, 0); ++j)
			{
				const size_t coID = m_rRefElem.id(d, i, 0, j);
				m_gloMid[d][i] += m_gloMid[0][coID];
			}

		// 	scale for correct averaging
			m_gloMid[d][i] *= 1./(m_rRefElem.num(d, i, 0));
		}
	}
	
	// set global integration points
	for (size_t i = 0; i < num_bf(); ++i)
	{
		// copy global corners of bf
		copy_global_corners(m_vBF[i]);
		
		if (dim == 1) // Edge
			{AveragePositions(m_vBF[i].globalIP, m_vBF[i].vGloPos, 2);}
		else if (dim == 2)	// Quadrilateral
			{AveragePositions(m_vBF[i].globalIP, m_vBF[i].vGloPos, 4);}
		else {UG_ASSERT(0, "Dimension higher than 2 not implemented.");}
	}
	
	// 	compute size of BFs
	for (size_t i = 0; i < num_bf(); ++i)
	{
	// 	copy global corners
		copy_global_corners(m_vBF[i]);

	// 	compute volume of bf
		m_vBF[i].vol = ElementSize<bf_type, worldDim>(m_vBF[i].vGloPos);
	}
	
	///////////////////////////
	// Copy ip pos in list
	///////////////////////////

	// 	loop Boundary Faces (BF)
	m_vGlobBFIP.clear();
	for (size_t i = 0; i < num_bf(); ++i)
	{
	//	get current BF
		const BF& rBF = bf(i);

	// 	loop ips
		for (size_t ip = 0; ip < rBF.num_ip(); ++ip)
		{
			m_vGlobBFIP.push_back(rBF.global_ip(ip));
		}
	}
}

//////////////////////
// FV1Geometry

template class FV1Geometry<Edge, 1>;
template class FV1Geometry<Edge, 2>;
template class FV1Geometry<Edge, 3>;

template class FV1Geometry<Triangle, 2>;
template class FV1Geometry<Triangle, 3>;

template class FV1Geometry<Quadrilateral, 2>;
template class FV1Geometry<Quadrilateral, 3>;

template class FV1Geometry<Tetrahedron, 3>;
template class FV1Geometry<Prism, 3>;
template class FV1Geometry<Pyramid, 3>;
template class FV1Geometry<Hexahedron, 3>;

//////////////////////
// DimFV1Geometry
template class DimFV1Geometry<1, 1>;
template class DimFV1Geometry<1, 2>;
template class DimFV1Geometry<1, 3>;

template class DimFV1Geometry<2, 2>;
template class DimFV1Geometry<2, 3>;

template class DimFV1Geometry<3, 3>;

//////////////////////
// Manifold
template class FV1ManifoldBoundary<Edge, 2>;
template class FV1ManifoldBoundary<Triangle, 3>;
template class FV1ManifoldBoundary<Quadrilateral, 3>;

} // end namespace ug

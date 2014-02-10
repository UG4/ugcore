/*
 * fvho_geom.cpp
 *
 *  Created on: 04.09.2010
 *      Author: andreasvogel
 */


#include "common/util/provider.h"
#include "fvho_geom.h"
#include "lib_disc/reference_element/reference_element.h"
#include "lib_disc/quadrature/quadrature_provider.h"
#include "lib_algebra/common/operations_vec.h"
#include <math.h>       /* pow */

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
	if (rRefElem.roid() == ROID_PYRAMID)
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

	if (rRefElem.roid() != ROID_PYRAMID)
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

	if (rRefElem.roid() != ROID_PYRAMID)
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

	if (rRefElem.roid() != ROID_PYRAMID || side != 0)
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

template <int dim>
void ComputeMultiIndicesOfSubElement(std::vector<MathVector<dim, int> >* vvMultiIndex,
                                     bool* vIsBndElem,
                                     std::vector<int>* vElemBndSide,
                                     std::vector<size_t>* vIndex,
                                     ReferenceObjectID roid,
                                     int p);

template <>
void ComputeMultiIndicesOfSubElement<1>(std::vector<MathVector<1, int> >* vvMultiIndex,
                                        bool* vIsBndElem,
                                        std::vector<int>* vElemBndSide,
                                        std::vector<size_t>* vIndex,
                                        ReferenceObjectID roid,
                                        int p)
{
//	switch for roid
	size_t se = 0;
	switch(roid)
	{
		case ROID_EDGE:
			for(int i = 0; i < p; ++i)
			{
				vvMultiIndex[se].resize(2);
				vvMultiIndex[se][0] = MathVector<1,int>(i);
				vvMultiIndex[se][1] = MathVector<1,int>(i+1);

				// reset bnd info
				vIsBndElem[se] = false;
				vElemBndSide[se].clear(); vElemBndSide[se].resize(3, -1);

				if(i==0)
				{
					vIsBndElem[se] = true;
					vElemBndSide[se][0] = 0;
				}
				if(i==p-1)
				{
					vIsBndElem[se] = true;
					vElemBndSide[se][1] = 1;
				}
				++se;
			}

			{
				FlexLagrangeLSFS<ReferenceEdge> set(p);
				for(size_t s = 0; s < se; ++s)
				{
					vIndex[s].resize(vvMultiIndex[s].size());
					for(size_t i = 0; i < vvMultiIndex[s].size(); ++i)
					{
						vIndex[s][i] = set.index(vvMultiIndex[s][i]);
					}
				}
			}
			break;
		default: UG_THROW("ReferenceElement not found.");
	}
}

template <>
void ComputeMultiIndicesOfSubElement<2>(std::vector<MathVector<2, int> >* vvMultiIndex,
                                        bool* vIsBndElem,
                                        std::vector<int>* vElemBndSide,
                                        std::vector<size_t>* vIndex,
                                        ReferenceObjectID roid,
                                        int p)
{
//	switch for roid
	size_t se = 0;
	switch(roid)
	{
		case ROID_TRIANGLE:
			for(int j = 0; j < p; ++j) { // y -direction
				for(int i = 0; i < p - j; ++i) { // x - direction
					vvMultiIndex[se].resize(3);
					vvMultiIndex[se][0] = MathVector<2,int>(i  , j  );
					vvMultiIndex[se][1] = MathVector<2,int>(i+1, j  );
					vvMultiIndex[se][2] = MathVector<2,int>(i  , j+1);

					// reset bnd info
					vIsBndElem[se] = false;
					vElemBndSide[se].clear(); vElemBndSide[se].resize(3, -1);

					if(i==0) // left
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][2] = 2;
					}
					if(j==0) // bottom
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][0] = 0;
					}
					if(i+j==p-1) // diag
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][1] = 1;
					}
					++se;
				}
			}

			for(int j = 1; j <= p; ++j) {
				for(int i = 1; i <= p - j; ++i) {
					vvMultiIndex[se].resize(3);
					vvMultiIndex[se][0] = MathVector<2,int>(i  , j  );
					vvMultiIndex[se][1] = MathVector<2,int>(i-1, j  );
					vvMultiIndex[se][2] = MathVector<2,int>(i  , j-1);

					// reset bnd info
					// all inner elems
					vIsBndElem[se] = false;
					vElemBndSide[se].clear(); vElemBndSide[se].resize(3, -1);
					++se;
				}
			}

			{
				FlexLagrangeLSFS<ReferenceTriangle> set(p);
				for(size_t s = 0; s < se; ++s)
				{
					vIndex[s].resize(vvMultiIndex[s].size());
					for(size_t i = 0; i < vvMultiIndex[s].size(); ++i)
					{
						vIndex[s][i] = set.index(vvMultiIndex[s][i]);
					}
				}
			}

			break;
		case ROID_QUADRILATERAL:
			for(int j = 0; j < p; ++j) {
				for(int i = 0; i < p; ++i) {
					vvMultiIndex[se].resize(4);
					vvMultiIndex[se][0] = MathVector<2,int>(i  , j  );
					vvMultiIndex[se][1] = MathVector<2,int>(i+1, j  );
					vvMultiIndex[se][2] = MathVector<2,int>(i+1, j+1);
					vvMultiIndex[se][3] = MathVector<2,int>(i  , j+1);

					// reset bnd info
					vIsBndElem[se] = false;
					vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

					if(i==0) // left
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][3] = 3;
					}
					if(i==p-1) // right
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][1] = 1;
					}
					if(j==0) // bottom
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][0] = 0;
					}
					if(j==p-1) // top
					{
						vIsBndElem[se] = true;
						vElemBndSide[se][2] = 2;
					}
					++se;
				}
			}

			{
				FlexLagrangeLSFS<ReferenceQuadrilateral> set(p);
				for(size_t s = 0; s < se; ++s)
				{
					vIndex[s].resize(vvMultiIndex[s].size());
					for(size_t i = 0; i < vvMultiIndex[s].size(); ++i)
					{
						vIndex[s][i] = set.index(vvMultiIndex[s][i]);
					}
				}
			}
			break;
		default: UG_THROW("ReferenceElement not found.");
	}

}

template <>
void ComputeMultiIndicesOfSubElement<3>(std::vector<MathVector<3, int> >* vvMultiIndex,
                                        bool* vIsBndElem,
                                        std::vector<int>* vElemBndSide,
                                        std::vector<size_t>* vIndex,
                                        ReferenceObjectID roid,
                                        int p)
{
//	switch for roid
	size_t se = 0;
	switch(roid)
	{
		case ROID_TETRAHEDRON:
			for(int k = 0; k < p; ++k) {
				for(int j = 0; j < p -k; ++j) {
					for(int i = 0; i < p -k -j; ++i) {
						vvMultiIndex[se].resize(4);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i+1, j  , k);
						vvMultiIndex[se][2] = MathVector<3,int>(i  , j+1, k);
						vvMultiIndex[se][3] = MathVector<3,int>(i  , j  , k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

						if(i==0) // left
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][2] = 2;
						}
						if(j==0) // front
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][3] = 3;
						}
						if(k==0) // bottom
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][0] = 0;
						}
						if(i+j+k==p-1) // diag
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][1] = 1;
						}
						++se;
					}
				}
			}
			//	build 4 tetrahedrons out of the remaining octogons
			for(int k = 0; k < p; ++k) {
				for(int j = 1; j < p -k; ++j) {
					for(int i = 0; i < p -k -j; ++i) {
						if(j >= 2){
							vvMultiIndex[se].resize(4);
							vvMultiIndex[se][0] = MathVector<3,int>(i+1, j-1, k);
							vvMultiIndex[se][1] = MathVector<3,int>(i+1, j-1, k+1);
							vvMultiIndex[se][2] = MathVector<3,int>(i  , j-1, k+1);
							vvMultiIndex[se][3] = MathVector<3,int>(i+1, j-2, k+1);

							// reset bnd info
							vIsBndElem[se] = false;
							vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

							++se;
						}

						vvMultiIndex[se].resize(4);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i  , j  , k+1);
						vvMultiIndex[se][2] = MathVector<3,int>(i  , j-1, k+1);
						vvMultiIndex[se][3] = MathVector<3,int>(i+1, j-1, k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

						if(i==0) // left
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][0] = 2;
						}
						++se;

						vvMultiIndex[se].resize(4);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i+1, j  , k);
						vvMultiIndex[se][2] = MathVector<3,int>(i+1, j-1, k+1);
						vvMultiIndex[se][3] = MathVector<3,int>(i+1, j-1, k);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

						if(k==0) // bottom
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][2] = 0;
						}
						++se;

						vvMultiIndex[se].resize(4);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i+1, j-1, k);
						vvMultiIndex[se][2] = MathVector<3,int>(i+1, j-1, k+1);
						vvMultiIndex[se][3] = MathVector<3,int>(i  , j-1, k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

						if(j==1) // front
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][1] = 3;
						}
						++se;

						vvMultiIndex[se].resize(4);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i+1, j-1, k+1);
						vvMultiIndex[se][2] = MathVector<3,int>(i+1, j  , k);
						vvMultiIndex[se][3] = MathVector<3,int>(i  , j  , k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(4, -1);

						if(i+j+k==p-1) // diag
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][1] = 1;
						}
						++se;
					}
				}
			}

			{
				FlexLagrangeLSFS<ReferenceTetrahedron> set(p);
				for(size_t s = 0; s < se; ++s)
				{
					vIndex[s].resize(vvMultiIndex[s].size());
					for(size_t i = 0; i < vvMultiIndex[s].size(); ++i)
					{
						vIndex[s][i] = set.index(vvMultiIndex[s][i]);
					}
				}
			}
			break;

		case ROID_PRISM:
			for(int k = 0; k < p; ++k) {
				for(int j = 0; j < p; ++j) {
					for(int i = 0; i < p - j; ++i) {
						vvMultiIndex[se].resize(6);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i+1, j  , k);
						vvMultiIndex[se][2] = MathVector<3,int>(i  , j+1, k);
						vvMultiIndex[se][3] = MathVector<3,int>(i  , j  , k+1);
						vvMultiIndex[se][4] = MathVector<3,int>(i+1, j  , k+1);
						vvMultiIndex[se][5] = MathVector<3,int>(i  , j+1, k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(5, -1);

						if(i==0) // left
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][3] = 3;
						}
						if(j==0) // front
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][1] = 1;
						}
						if(i+j==p-1) // diag
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][2] = 2;
						}
						if(k==0) // bottom
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][0] = 0;
						}
						if(k==p-1) // top
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][4] = 4;
						}
						++se;
					}
				}
			}
			for(int k = 0; k < p; ++k) {
				for(int j = 1; j <= p; ++j) {
					for(int i = 1; i <= p - j; ++i) {
						vvMultiIndex[se].resize(6);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  ,k);
						vvMultiIndex[se][1] = MathVector<3,int>(i-1, j  ,k);
						vvMultiIndex[se][2] = MathVector<3,int>(i  , j-1,k);
						vvMultiIndex[se][3] = MathVector<3,int>(i  , j  ,k+1);
						vvMultiIndex[se][4] = MathVector<3,int>(i-1, j  ,k+1);
						vvMultiIndex[se][5] = MathVector<3,int>(i  , j-1,k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(5, -1);

						if(k==0) // bottom
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][0] = 0;
						}
						if(k==p-1) // top
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][4] = 4;
						}
						++se;
					}
				}
			}

			{
				FlexLagrangeLSFS<ReferencePrism> set(p);
				for(size_t s = 0; s < se; ++s)
				{
					vIndex[s].resize(vvMultiIndex[s].size());
					for(size_t i = 0; i < vvMultiIndex[s].size(); ++i)
					{
						vIndex[s][i] = set.index(vvMultiIndex[s][i]);
					}
				}
			}

			break;
		case ROID_HEXAHEDRON:
			for(int k = 0; k < p; ++k) {
				for(int j = 0; j < p; ++j) {
					for(int i = 0; i < p; ++i) {
						vvMultiIndex[se].resize(8);
						vvMultiIndex[se][0] = MathVector<3,int>(i  , j  , k);
						vvMultiIndex[se][1] = MathVector<3,int>(i+1, j  , k);
						vvMultiIndex[se][2] = MathVector<3,int>(i+1, j+1, k);
						vvMultiIndex[se][3] = MathVector<3,int>(i  , j+1, k);
						vvMultiIndex[se][4] = MathVector<3,int>(i  , j  , k+1);
						vvMultiIndex[se][5] = MathVector<3,int>(i+1, j  , k+1);
						vvMultiIndex[se][6] = MathVector<3,int>(i+1, j+1, k+1);
						vvMultiIndex[se][7] = MathVector<3,int>(i  , j+1, k+1);

						// reset bnd info
						vIsBndElem[se] = false;
						vElemBndSide[se].clear(); vElemBndSide[se].resize(6, -1);

						if(i==0) // left
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][4] = 4;
						}
						if(i==p-1) //right
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][2] = 2;
						}
						if(j==0) // front
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][1] = 1;
						}
						if(j==p-1) // back
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][3] = 3;
						}
						if(k==0) // bottom
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][0] = 0;
						}
						 if(k==p-1) // top
						{
							vIsBndElem[se] = true;
							vElemBndSide[se][5] = 5;
						}
						++se;
					}
				}
			}

			{
				FlexLagrangeLSFS<ReferenceHexahedron> set(p);
				for(size_t s = 0; s < se; ++s)
				{
					vIndex[s].resize(vvMultiIndex[s].size());
					for(size_t i = 0; i < vvMultiIndex[s].size(); ++i)
					{
						vIndex[s][i] = set.index(vvMultiIndex[s][i]);
					}
				}
			}

			break;
		default: UG_THROW("ReferenceElement not found.");
	}

}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// FV Geometry for Reference Element Type (all order, FVHO)
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int TOrder, typename TElem, int TWorldDim, int TQuadOrder>
FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>::
FVGeometry()
	: m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()),
	  m_rTrialSpace(Provider<local_shape_fct_set_type>::get()),
	  m_rSCVFQuadRule(Provider<scvf_quad_rule_type>::get()),
	  m_rSCVQuadRule(Provider<scv_quad_rule_type>::get())
{
	update_local_data();
}

template <int TOrder, typename TElem, int TWorldDim, int TQuadOrder>
void FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>::
update_local(ReferenceObjectID roid, const LFEID& lfeID,
                  size_t quadOrder)
{
	if(roid != geometry_traits<TElem>::REFERENCE_OBJECT_ID)
	{
		UG_THROW("FVGeometry::update: Geometry only for "
				<<geometry_traits<TElem>::REFERENCE_OBJECT_ID<<", but "
				<<roid<<" requested.");
	}
	if(lfeID.type() != LFEID::LAGRANGE)
	{
		UG_THROW("FVGeometry::update: Geometry only for shape type"
				<<"Lagrange"<<", but "<<lfeID.type()<<" requested.");
	}
	if(lfeID.order() != TOrder)
	{
		UG_THROW("FVGeometry::update: Geometry only for shape order"
				<<TOrder<<", but "<<lfeID.order()<<" requested.");
	}
	if(quadOrder > TQuadOrder)
	{
		UG_THROW("FVGeometry::update: Geometry only for scvf integration order "
				<< TQuadOrder<<", but order "<<quadOrder<<" requested.");
	}
}


template <int TOrder, typename TElem, int TWorldDim, int TQuadOrder>
void FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>::
update_local_data()
{
//	get reference object id
	ReferenceObjectID roid = ref_elem_type::REFERENCE_OBJECT_ID;

//	determine corners of sub elements
	bool vIsBndElem[numSubElem];
	std::vector<int> vElemBndSide[numSubElem];
	std::vector<size_t> vIndex[numSubElem];
	std::vector<MathVector<dim,int> > vMultiIndex[numSubElem];

	ComputeMultiIndicesOfSubElement<dim>(vMultiIndex, vIsBndElem,
										 vElemBndSide, vIndex, roid, p);

//	directions of counting
	MathVector<dim> direction[dim];
	for(int i = 0; i < dim; ++i){direction[i] = 0.0; direction[i][i] = 1.0;}

	for(size_t se = 0; se < numSubElem; ++se)
	{
		for(int co = 0; co < ref_elem_type::numCorners; ++co)
		{
		//	compute corners of sub elem in local coordinates
			MathVector<dim> pos; pos = 0.0;
			for(int i = 0; i < dim; ++i)
			{
				const number frac = vMultiIndex[se][co][i] / ((number)p);
				VecScaleAppend(pos, frac, direction[i]);
			}
			m_vSubElem[se].vvLocMid[0][co] = pos;

		//	get multi index for corner
			m_vSubElem[se].vDoFID[co] = vIndex[se][co];
		}

	//	remember if boundary element
		m_vSubElem[se].isBndElem = vIsBndElem[se];

	//	remember boundary sides
		m_vSubElem[se].vElemBndSide = vElemBndSide[se];
	}

//	compute mid points for all sub elements
	for(size_t se = 0; se < numSubElem; ++se)
		ComputeMidPoints<dim, ref_elem_type, maxMid>
				(m_rRefElem, m_vSubElem[se].vvLocMid[0], m_vSubElem[se].vvLocMid);

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one edge of the sub-element
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	get corresponding subelement
		const size_t se = i / numSCVFPerSubElem;
		const size_t locSCVF = i % numSCVFPerSubElem;

	//	this scvf separates the given nodes
		const size_t locFrom =  m_rRefElem.id(1, locSCVF, 0, 0);
		const size_t locTo =  m_rRefElem.id(1, locSCVF, 0, 1);

		m_vSCVF[i].From = m_vSubElem[se].vDoFID[locFrom];
		m_vSCVF[i].To = m_vSubElem[se].vDoFID[locTo];

	//	compute mid ids of the scvf
		ComputeSCVFMidID(m_rRefElem, m_vSCVF[i].vMidID, locSCVF);

	// 	copy local corners of scvf
		CopyCornerByMidID<dim, maxMid>
			(m_vSCVF[i].vLocPos, m_vSCVF[i].vMidID, m_vSubElem[se].vvLocMid, SCVF::numCo);

	// 	compute integration points
		m_vSCVF[i].vWeight = m_rSCVFQuadRule.weights();
		ReferenceMapping<scvf_type, dim> map(m_vSCVF[i].vLocPos);
		for(size_t ip = 0; ip < m_rSCVFQuadRule.size(); ++ip)
			map.local_to_global(m_vSCVF[i].vLocalIP[ip], m_rSCVFQuadRule.point(ip));
	}


// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the sub-element
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	get corresponding subelement
		const size_t se = i / numSCVPerSubElem;
		const size_t locSCV = i % numSCVPerSubElem;

	//	store associated node
		m_vSCV[i].nodeId = m_vSubElem[se].vDoFID[locSCV];;

	//	compute mid ids scv
		ComputeSCVMidID(m_rRefElem, m_vSCV[i].vMidID, locSCV);

	// 	copy local corners of scv
		CopyCornerByMidID<dim, maxMid>
			(m_vSCV[i].vLocPos, m_vSCV[i].vMidID, m_vSubElem[se].vvLocMid, m_vSCV[i].num_corners());

	// 	compute integration points
		m_vSCV[i].vWeight = m_rSCVQuadRule.weights();
		ReferenceMapping<scv_type, dim> map(m_vSCV[i].vLocPos);
		for(size_t ip = 0; ip < m_rSCVQuadRule.size(); ++ip)
			map.local_to_global(m_vSCV[i].vLocalIP[ip], m_rSCVQuadRule.point(ip));
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < m_vSCVF[i].num_ip(); ++ip)
		{
			m_rTrialSpace.shapes(&(m_vSCVF[i].vvShape[ip][0]), m_vSCVF[i].local_ip(ip));
			m_rTrialSpace.grads(&(m_vSCVF[i].vvLocalGrad[ip][0]), m_vSCVF[i].local_ip(ip));
		}

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < m_vSCV[i].num_ip(); ++ip)
		{
			m_rTrialSpace.shapes(&(m_vSCV[i].vvShape[ip][0]), m_vSCV[i].local_ip(ip));
			m_rTrialSpace.grads(&(m_vSCV[i].vvLocalGrad[ip][0]), m_vSCV[i].local_ip(ip));
		}

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	size_t allIP = 0;
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < m_vSCVF[i].num_ip(); ++ip)
			m_vLocSCVF_IP[allIP++] = scvf(i).local_ip(ip);

// 	copy ip positions in a list for Sub Control Volumes (SCV)
	allIP = 0;
	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < m_vSCV[i].num_ip(); ++ip)
			m_vLocSCV_IP[allIP++] = scv(i).local_ip(ip);
}


/// update data for given element
template <int TOrder, typename TElem, int TWorldDim, int TQuadOrder>
void FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>::
update(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

//	update reference mapping
	m_rMapping.update(vCornerCoords);

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	map local corners of scvf to global
		for(size_t co = 0; co < m_vSCVF[i].num_corners(); ++co)
			m_rMapping.local_to_global(m_vSCVF[i].vGloPos[co], m_vSCVF[i].vLocPos[co]);

	//	map local ips of scvf to global
		for(size_t ip = 0; ip < m_vSCVF[i].num_ip(); ++ip)
			m_rMapping.local_to_global(m_vSCVF[i].vGlobalIP[ip], m_vSCVF[i].local_ip(ip));

	// 	normal on scvf
		traits::NormalOnSCVF(m_vSCVF[i].Normal, m_vSCVF[i].vGloPos, vCornerCoords);
		VecNormalize(m_vSCVF[i].Normal, m_vSCVF[i].Normal);

		ReferenceMapping<scvf_type, dim> map(m_vSCVF[i].vGloPos);
		for(size_t ip = 0; ip < m_rSCVFQuadRule.size(); ++ip)
			m_vSCVF[i].vDetJMap[ip] = map.sqrt_gram_det(m_rSCVFQuadRule.point(ip));
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	map local corners of scvf to global
		for(size_t co = 0; co < m_vSCV[i].num_corners(); ++co)
			m_rMapping.local_to_global(m_vSCV[i].vGloPos[co], m_vSCV[i].vLocPos[co]);

	//	map local ips of scvf to global
		for(size_t ip = 0; ip < m_vSCV[i].num_ip(); ++ip)
			m_rMapping.local_to_global(m_vSCV[i].vGlobalIP[ip], m_vSCV[i].local_ip(ip));
	}

//	if mapping is linear, compute jacobian only once and copy
	if(ReferenceMapping<ref_elem_type, worldDim>::isLinear)
	{
		MathMatrix<worldDim,dim> JtInv;
		m_rMapping.jacobian_transposed_inverse(JtInv, m_vSCVF[0].local_ip(0));
		const number detJ = m_rMapping.sqrt_gram_det(m_vSCVF[0].local_ip(0));
		for(size_t i = 0; i < num_scvf(); ++i)
			for(size_t ip = 0; ip < scvf(i).num_ip(); ++ip)
			{
				m_vSCVF[i].vJtInv[ip] = JtInv;
				m_vSCVF[i].vDetJ[ip] = detJ;
			}

		for(size_t i = 0; i < num_scv(); ++i)
			for(size_t ip = 0; ip < scv(i).num_ip(); ++ip)
			{
				m_vSCV[i].vJtInv[ip] = JtInv;
				m_vSCV[i].vDetJ[ip] = detJ;
			}
	}
//	else compute jacobian for each integration point
	else
	{
		for(size_t i = 0; i < num_scvf(); ++i)
			for(size_t ip = 0; ip < m_vSCVF[i].num_ip(); ++ip)
			{
				m_rMapping.jacobian_transposed_inverse(m_vSCVF[i].vJtInv[ip], m_vSCVF[i].local_ip(ip));
				m_vSCVF[i].vDetJ[ip] = m_rMapping.sqrt_gram_det(m_vSCVF[i].local_ip(ip));
			}

		for(size_t i = 0; i < num_scv(); ++i)
			for(size_t ip = 0; ip < m_vSCV[i].num_ip(); ++ip)
			{
				m_rMapping.jacobian_transposed_inverse(m_vSCV[i].vJtInv[ip], m_vSCV[i].local_ip(ip));
				m_vSCV[i].vDetJ[ip] = m_rMapping.sqrt_gram_det(m_vSCV[i].local_ip(ip));
			}
	}

	for(size_t i = 0; i < num_scv(); ++i)
	{
		ReferenceMapping<scv_type, dim> map(m_vSCV[i].vGloPos);
		for(size_t ip = 0; ip < m_rSCVQuadRule.size(); ++ip)
			m_vSCV[i].vDetJMap[ip] = map.sqrt_gram_det(m_rSCVQuadRule.point(ip));
	}

//	compute global gradients
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < scvf(i).num_ip(); ++ip)
			for(size_t sh = 0 ; sh < nsh; ++sh)
				MatVecMult(m_vSCVF[i].vvGlobalGrad[ip][sh], m_vSCVF[i].vJtInv[ip], m_vSCVF[i].vvLocalGrad[ip][sh]);

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < scv(i).num_ip(); ++ip)
			for(size_t sh = 0 ; sh < nsh; ++sh)
				MatVecMult(m_vSCV[i].vvGlobalGrad[ip][sh], m_vSCV[i].vJtInv[ip], m_vSCV[i].vvLocalGrad[ip][sh]);

// 	Copy ip pos in list for SCVF
	size_t allIP = 0;
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < scvf(i).num_ip(); ++ip)
			m_vGlobSCVF_IP[allIP++] = scvf(i).global_ip(ip);

	allIP = 0;
	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < scv(i).num_ip(); ++ip)
			m_vGlobSCV_IP[allIP++] = scv(i).global_ip(ip);

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}

template <int TOrder, typename TElem, int TWorldDim, int TQuadOrder>
void FVGeometry<TOrder, TElem, TWorldDim, TQuadOrder>::
update_boundary_faces(GeometricObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

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

//	update reference mapping
	m_rMapping.update(vCornerCoords);

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

	//	loop subelements
		for(size_t se = 0; se < numSubElem; ++se)
		{
		//	skip inner sub elements
			if(!m_vSubElem[se].isBndElem) continue;

		//	loop sides of element
			for(size_t side = 0; side < m_vSubElem[se].vElemBndSide.size(); ++side)
			{
			//	get whole element bnd side
				const int elemBndSide = m_vSubElem[se].vElemBndSide[side];

			//	skip non boundary sides
				if(elemBndSide == -1 || vSubsetIndex[elemBndSide] != bndIndex) continue;

			//	number of corners of side
				const int coOfSide = m_rRefElem.num(dim-1, elemBndSide, 0);

			//	resize vector
				vBF.resize(curr_bf + coOfSide);

			//	loop corners
				for(int co = 0; co < coOfSide; ++co)
				{
				//	get current bf
					BF& bf = vBF[curr_bf];

				//	set node id == scv this bf belongs to
					const int refNodeId = m_rRefElem.id(dim-1, elemBndSide, 0, co);
					bf.nodeId = m_vSubElem[se].vDoFID[refNodeId];

				//	Compute MidID for BF
					ComputeBFMidID(m_rRefElem, elemBndSide, bf.vMidID, co);

				// 	copy corners of bf
					CopyCornerByMidID<dim, maxMid>
						(bf.vLocPos, bf.vMidID, m_vSubElem[se].vvLocMid, BF::numCo);
//					CopyCornerByMidID<worldDim, maxMid>
//						(bf.vGloPos, bf.vMidID, m_vSubElem[se].vvGloMid, BF::numCo);
				// 	compute global corners of bf
					for(size_t i = 0; i < bf.num_corners(); ++i)
						m_rMapping.local_to_global(bf.vGloPos[i], bf.vLocPos[i]);


				// 	normal on scvf
					traits::NormalOnSCVF(bf.Normal, bf.vGloPos, m_vSubElem[se].vvGloMid[0]);

				//	compute local integration points
					bf.vWeight = m_rSCVFQuadRule.weights();
					ReferenceMapping<scvf_type, dim> map(bf.vLocPos);
					for(size_t ip = 0; ip < m_rSCVFQuadRule.size(); ++ip)
						map.local_to_global(bf.vLocalIP[ip], m_rSCVFQuadRule.point(ip));

				//	compute global integration points
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
						m_rMapping.local_to_global(bf.vGlobalIP[ip], bf.vLocalIP[ip]);

				//	compute volume
					bf.Vol = VecTwoNorm(bf.Normal);

				//	compute shapes and gradients
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
					{
						m_rTrialSpace.shapes(&(bf.vvShape[ip][0]), bf.local_ip(ip));
						m_rTrialSpace.grads(&(bf.vvLocalGrad[ip][0]), bf.local_ip(ip));

						m_rMapping.jacobian_transposed_inverse(bf.vJtInv[ip], bf.local_ip(ip));
						bf.vDetJ[ip] = m_rMapping.sqrt_gram_det(bf.local_ip(ip));
					}

				//	compute global gradient
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
						for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
							MatVecMult(bf.vvGlobalGrad[ip][sh],
							           bf.vJtInv[ip], bf.vvLocalGrad[ip][sh]);

				//	increase curr_bf
					++curr_bf;
				}
			}
		}
	}
}


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// FV Geometry (all order, FVHO)   DIM FV
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <int TDim, int TWorldDim>
DimFVGeometry<TDim, TWorldDim>::
DimFVGeometry()	: m_pElem(NULL) {}



template <int TDim, int TWorldDim>
void DimFVGeometry<TDim, TWorldDim>::
update_local(ReferenceObjectID roid, const LFEID& lfeID, size_t orderQuad)
{
//	save setting we prepare the local data for
	m_roid = roid;
	m_lfeID = lfeID;
	m_orderShape = m_lfeID.order();
	m_quadOrderSCVF = orderQuad;
	m_quadOrderSCV = orderQuad;

//	resize sub elements
	m_numSubElem = (size_t) std::pow((double) m_orderShape, dim);
	m_vSubElem.resize(m_numSubElem);

//	get the multi indices for the sub elements and the boundary flags
	bool* vIsBndElem = new bool[m_numSubElem];
	std::vector<std::vector<int> > vElemBndSide(m_numSubElem);
	std::vector<std::vector<MathVector<dim,int> > > vMultiIndex(m_numSubElem);
	std::vector<std::vector<size_t> > vIndex(m_numSubElem);
	ComputeMultiIndicesOfSubElement<dim>(&vMultiIndex[0], &vIsBndElem[0],
	                                     &vElemBndSide[0], &vIndex[0], m_roid, m_orderShape);

//	get reference element
	try{
	const DimReferenceElement<dim>& rRefElem
				= ReferenceElementProvider::get<dim>(m_roid);

//	directions of counting
	MathVector<dim> direction[dim];
	for(int i = 0; i < dim; ++i){direction[i] = 0.0; direction[i][i] = 1.0;}

	for(size_t se = 0; se < m_numSubElem; ++se)
	{
		for(size_t co = 0; co < vMultiIndex[se].size(); ++co)
		{
		//	compute corners of sub elem in local coordinates
			MathVector<dim> pos; pos = 0.0;
			for(int i = 0; i < dim; ++i)
			{
				const number frac = vMultiIndex[se][co][i] / ((number)m_orderShape);
				VecScaleAppend(pos, frac, direction[i]);
			}
			m_vSubElem[se].vvLocMid[0][co] = pos;

		//	get multi index for corner
			m_vSubElem[se].vDoFID[co] = vIndex[se][co];
		}

	//	remember if boundary element
		m_vSubElem[se].isBndElem = vIsBndElem[se];

	//	remember boundary sides
		m_vSubElem[se].vElemBndSide = vElemBndSide[se];
	}

	delete[] vIsBndElem;

//	compute mid points for all sub elements
	for(size_t se = 0; se < m_numSubElem; ++se)
		ComputeMidPoints<dim, DimReferenceElement<dim>, maxMid>
				(rRefElem, m_vSubElem[se].vvLocMid[0], m_vSubElem[se].vvLocMid);

//	number of scvf/scv per subelem
	m_numSCVFPerSubElem = rRefElem.num(1);
	m_numSCVPerSubElem = rRefElem.num(0);

	m_numSCVF = m_numSCVFPerSubElem * m_numSubElem;
	m_numSCV = m_numSCVPerSubElem * m_numSubElem;

	m_vSCVF.resize(m_numSCVF);
	m_vSCV.resize(m_numSCV);


//	get trial space
	const LocalShapeFunctionSet<dim>& rTrialSpace =
		LocalFiniteElementProvider::get<dim>(m_roid, m_lfeID);

//	request for quadrature rule
	const ReferenceObjectID scvfRoid = scvf_type::REFERENCE_OBJECT_ID;
	const QuadratureRule<dim-1>& rSCVFQuadRule
			= QuadratureRuleProvider<dim-1>::get(scvfRoid, m_quadOrderSCVF);

	const int nipSCVF = rSCVFQuadRule.size();
	m_numSCVFIP = m_numSCVF * nipSCVF;

// 	set up local informations for SubControlVolumeFaces (scvf)
// 	each scvf is associated to one edge of the sub-element
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	get corresponding subelement
		const size_t se = i / m_numSCVFPerSubElem;
		const size_t locSCVF = i % m_numSCVFPerSubElem;

	//	this scvf separates the given nodes
		const size_t locFrom =  rRefElem.id(1, locSCVF, 0, 0);
		const size_t locTo =  rRefElem.id(1, locSCVF, 0, 1);

		m_vSCVF[i].From = m_vSubElem[se].vDoFID[locFrom];
		m_vSCVF[i].To = m_vSubElem[se].vDoFID[locTo];

	//	compute mid ids of the scvf
		ComputeSCVFMidID(rRefElem, m_vSCVF[i].vMidID, locSCVF);

	// 	copy local corners of scvf
		CopyCornerByMidID<dim, maxMid>
			(m_vSCVF[i].vLocPos, m_vSCVF[i].vMidID, m_vSubElem[se].vvLocMid, SCVF::numCo);

	// 	compute integration points
		m_vSCVF[i].vWeight = rSCVFQuadRule.weights();
		m_vSCVF[i].nip = nipSCVF;

		m_vSCVF[i].vLocalIP.resize(nipSCVF);
		m_vSCVF[i].vGlobalIP.resize(nipSCVF);

		m_vSCVF[i].vvShape.resize(nipSCVF);
		m_vSCVF[i].vvLocalGrad.resize(nipSCVF);
		m_vSCVF[i].vvGlobalGrad.resize(nipSCVF);
		m_vSCVF[i].vJtInv.resize(nipSCVF);
		m_vSCVF[i].vDetJ.resize(nipSCVF);
		m_vSCVF[i].vDetJMap.resize(nipSCVF);

		m_vSCVF[i].nsh = rTrialSpace.num_sh();

		ReferenceMapping<scvf_type, dim> map(m_vSCVF[i].vLocPos);
		for(size_t ip = 0; ip < rSCVFQuadRule.size(); ++ip)
			map.local_to_global(m_vSCVF[i].vLocalIP[ip], rSCVFQuadRule.point(ip));
	}


//	request for quadrature rule
	static const ReferenceObjectID scvRoid = scv_type::REFERENCE_OBJECT_ID;
	const QuadratureRule<dim>& rSCVQuadRule
			= QuadratureRuleProvider<dim>::get(scvRoid, m_quadOrderSCV);

	const int nipSCV = rSCVQuadRule.size();
	m_numSCVIP = m_numSCV * nipSCV;

// 	set up local informations for SubControlVolumes (scv)
// 	each scv is associated to one corner of the sub-element
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	get corresponding subelement
		const size_t se = i / m_numSCVPerSubElem;
		const size_t locSCV = i % m_numSCVPerSubElem;

	//	store associated node
		m_vSCV[i].nodeId = m_vSubElem[se].vDoFID[locSCV];;

	//	compute mid ids scv
		ComputeSCVMidID(rRefElem, m_vSCV[i].vMidID, locSCV);

	// 	copy local corners of scv
		CopyCornerByMidID<dim, maxMid>
			(m_vSCV[i].vLocPos, m_vSCV[i].vMidID, m_vSubElem[se].vvLocMid, m_vSCV[i].num_corners());

	// 	compute integration points
		m_vSCV[i].vWeight = rSCVQuadRule.weights();
		m_vSCV[i].nip = nipSCV;

		m_vSCV[i].vLocalIP.resize(nipSCV);
		m_vSCV[i].vGlobalIP.resize(nipSCV);

		m_vSCV[i].vvShape.resize(nipSCV);
		m_vSCV[i].vvLocalGrad.resize(nipSCV);
		m_vSCV[i].vvGlobalGrad.resize(nipSCV);
		m_vSCV[i].vJtInv.resize(nipSCV);
		m_vSCV[i].vDetJ.resize(nipSCV);
		m_vSCV[i].vDetJMap.resize(nipSCV);

		m_vSCV[i].nsh = rTrialSpace.num_sh();

		ReferenceMapping<scv_type, dim> map(m_vSCV[i].vLocPos);
		for(size_t ip = 0; ip < rSCVQuadRule.size(); ++ip)
			map.local_to_global(m_vSCV[i].vLocalIP[ip], rSCVQuadRule.point(ip));
	}

	/////////////////////////
	// Shapes and Derivatives
	/////////////////////////

	m_nsh = rTrialSpace.num_sh();

	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < m_vSCVF[i].num_ip(); ++ip)
		{
			m_vSCVF[i].vvShape[ip].resize(m_vSCVF[i].nsh);
			m_vSCVF[i].vvLocalGrad[ip].resize(m_vSCVF[i].nsh);
			m_vSCVF[i].vvGlobalGrad[ip].resize(m_vSCVF[i].nsh);

			rTrialSpace.shapes(&(m_vSCVF[i].vvShape[ip][0]), m_vSCVF[i].local_ip(ip));
			rTrialSpace.grads(&(m_vSCVF[i].vvLocalGrad[ip][0]), m_vSCVF[i].local_ip(ip));
		}

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < m_vSCV[i].num_ip(); ++ip)
		{
			m_vSCV[i].vvShape[ip].resize(m_vSCV[i].nsh);
			m_vSCV[i].vvLocalGrad[ip].resize(m_vSCV[i].nsh);
			m_vSCV[i].vvGlobalGrad[ip].resize(m_vSCV[i].nsh);

			rTrialSpace.shapes(&(m_vSCV[i].vvShape[ip][0]), m_vSCV[i].local_ip(ip));
			rTrialSpace.grads(&(m_vSCV[i].vvLocalGrad[ip][0]), m_vSCV[i].local_ip(ip));
		}

	}
	UG_CATCH_THROW("DimFV1Geometry: update failed.");

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	m_vLocSCVF_IP.resize(m_numSCVFIP);
	m_vGlobSCVF_IP.resize(m_numSCVFIP);
	size_t allIP = 0;
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < m_vSCVF[i].num_ip(); ++ip)
			m_vLocSCVF_IP[allIP++] = scvf(i).local_ip(ip);

// 	copy ip positions in a list for Sub Control Volumes (SCV)
	m_vLocSCV_IP.resize(m_numSCVIP);
	m_vGlobSCV_IP.resize(m_numSCVIP);
	allIP = 0;
	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < m_vSCV[i].num_ip(); ++ip)
			m_vLocSCV_IP[allIP++] = scv(i).local_ip(ip);
}


/// update data for given element
template <int TDim, int TWorldDim>
void DimFVGeometry<TDim, TWorldDim>::
update(GeometricObject* pElem, const MathVector<worldDim>* vCornerCoords,
       const LFEID& lfeID, size_t quadOrder,
       const ISubsetHandler* ish)
{
// 	If already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

//	get reference element type
	ReferenceObjectID roid = (ReferenceObjectID)pElem->reference_object_id();

//	if already prepared for this roid, skip update of local values
	if(m_roid != roid || lfeID != m_lfeID ||
	   (int)quadOrder != m_quadOrderSCVF || (int)quadOrder != m_quadOrderSCV)
			update_local(roid, lfeID, quadOrder);

//	get reference element mapping
	try{
	DimReferenceMapping<dim, worldDim>& rMapping
		= ReferenceMappingProvider::get<dim, worldDim>(roid);

//	update reference mapping
	rMapping.update(vCornerCoords);

// 	compute global informations for scvf
	for(size_t i = 0; i < num_scvf(); ++i)
	{
	//	map local corners of scvf to global
		for(size_t co = 0; co < m_vSCVF[i].num_corners(); ++co)
			rMapping.local_to_global(m_vSCVF[i].vGloPos[co], m_vSCVF[i].vLocPos[co]);

	// 	normal on scvf
		traits::NormalOnSCVF(m_vSCVF[i].Normal, m_vSCVF[i].vGloPos, vCornerCoords);
		VecNormalize(m_vSCVF[i].Normal, m_vSCVF[i].Normal);

		static const ReferenceObjectID scvfRoid = scvf_type::REFERENCE_OBJECT_ID;
		const QuadratureRule<dim-1>& rSCVFQuadRule
				= QuadratureRuleProvider<dim-1>::get(scvfRoid, m_quadOrderSCVF);
		ReferenceMapping<scvf_type, worldDim> map(m_vSCVF[i].vGloPos);

		for(size_t ip = 0; ip < rSCVFQuadRule.size(); ++ip){
			m_vSCVF[i].vDetJMap[ip] = map.sqrt_gram_det(rSCVFQuadRule.point(ip));
			map.local_to_global(m_vSCVF[i].vGlobalIP[ip], rSCVFQuadRule.point(ip));
		}

		rMapping.jacobian_transposed_inverse(&m_vSCVF[i].vJtInv[0], &m_vSCVF[i].vLocalIP[0], m_vSCVF[i].num_ip());
		rMapping.sqrt_gram_det(&m_vSCVF[i].vDetJ[0], &m_vSCVF[i].vLocalIP[0], m_vSCVF[i].num_ip());
	}

// 	compute size of scv
	for(size_t i = 0; i < num_scv(); ++i)
	{
	//	map local corners of scvf to global
		rMapping.local_to_global(&m_vSCV[i].vGloPos[0], &m_vSCV[i].vLocPos[0], m_vSCV[i].num_corners());

		rMapping.jacobian_transposed_inverse(&m_vSCV[i].vJtInv[0], &m_vSCV[i].vLocalIP[0], m_vSCV[i].num_ip());
		rMapping.sqrt_gram_det(&m_vSCV[i].vDetJ[0], &m_vSCV[i].vLocalIP[0], m_vSCV[i].num_ip());


		static const ReferenceObjectID scvRoid = scv_type::REFERENCE_OBJECT_ID;
		const QuadratureRule<dim>& rSCVQuadRule
				= QuadratureRuleProvider<dim>::get(scvRoid, m_quadOrderSCV);

		ReferenceMapping<scv_type, worldDim> map(m_vSCV[i].vGloPos);
		for(size_t ip = 0; ip < rSCVQuadRule.size(); ++ip){
			m_vSCV[i].vDetJMap[ip] = map.sqrt_gram_det(rSCVQuadRule.point(ip));
			map.local_to_global(m_vSCV[i].vGlobalIP[ip], rSCVQuadRule.point(ip));
		}
	}

	}
	UG_CATCH_THROW("DimFVGeometry: update failed.");

//	compute global gradients
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < scvf(i).num_ip(); ++ip)
			for(size_t sh = 0 ; sh < m_vSCVF[i].nsh; ++sh)
				MatVecMult(m_vSCVF[i].vvGlobalGrad[ip][sh], m_vSCVF[i].vJtInv[ip], m_vSCVF[i].vvLocalGrad[ip][sh]);

	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < scv(i).num_ip(); ++ip)
			for(size_t sh = 0 ; sh < m_vSCV[i].nsh; ++sh)
				MatVecMult(m_vSCV[i].vvGlobalGrad[ip][sh], m_vSCV[i].vJtInv[ip], m_vSCV[i].vvLocalGrad[ip][sh]);

// 	Copy ip pos in list for SCVF
	size_t allIP = 0;
	for(size_t i = 0; i < num_scvf(); ++i)
		for(size_t ip = 0; ip < scvf(i).num_ip(); ++ip)
			m_vGlobSCVF_IP[allIP++] = scvf(i).global_ip(ip);

	allIP = 0;
	for(size_t i = 0; i < num_scv(); ++i)
		for(size_t ip = 0; ip < scv(i).num_ip(); ++ip)
			m_vGlobSCV_IP[allIP++] = scv(i).global_ip(ip);

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}

template <int TDim, int TWorldDim>
void DimFVGeometry<TDim, TWorldDim>::
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

//	get reference element mapping
	try{
	DimReferenceMapping<dim, worldDim>& rMapping
		= ReferenceMappingProvider::get<dim, worldDim>(m_roid);

	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

	const ReferenceObjectID scvfRoid = scvf_type::REFERENCE_OBJECT_ID;
	const QuadratureRule<dim-1>& rSCVFQuadRule
			= QuadratureRuleProvider<dim-1>::get(scvfRoid, m_quadOrderSCVF);

	const LocalShapeFunctionSet<dim>& rTrialSpace =
		LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, dim, m_orderShape));

//	update reference mapping
	rMapping.update(vCornerCoords);

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

	//	loop subelements
		for(size_t se = 0; se < m_numSubElem; ++se)
		{
		//	skip inner sub elements
			if(!m_vSubElem[se].isBndElem) continue;

		//	loop sides of element
			for(size_t side = 0; side < m_vSubElem[se].vElemBndSide.size(); ++side)
			{
			//	get whole element bnd side
				const int elemBndSide = m_vSubElem[se].vElemBndSide[side];

			//	skip non boundary sides
				if(elemBndSide == -1 || vSubsetIndex[elemBndSide] != bndIndex) continue;

			//	number of corners of side
				const int coOfSide = rRefElem.num(dim-1, elemBndSide, 0);

			//	resize vector
				vBF.resize(curr_bf + coOfSide);

			//	loop corners
				for(int co = 0; co < coOfSide; ++co)
				{
				//	get current bf
					BF& bf = vBF[curr_bf];

				//	set node id == scv this bf belongs to
					const int refNodeId = rRefElem.id(dim-1, elemBndSide, 0, co);
					bf.nodeId = m_vSubElem[se].vDoFID[refNodeId];

				//	Compute MidID for BF
					ComputeBFMidID(rRefElem, elemBndSide, bf.vMidID, co);

				// 	copy corners of bf
					CopyCornerByMidID<dim, maxMid>
						(bf.vLocPos, bf.vMidID, m_vSubElem[se].vvLocMid, BF::numCo);
//					CopyCornerByMidID<worldDim, maxMid>
//						(bf.vGloPos, bf.vMidID, m_vSubElem[se].vvGloMid, BF::numCo);				// 	compute global corners of bf
					for(size_t i = 0; i < bf.num_corners(); ++i)
						rMapping.local_to_global(bf.vGloPos[i], bf.vLocPos[i]);

				// 	normal on scvf
					traits::NormalOnSCVF(bf.Normal, bf.vGloPos, m_vSubElem[se].vvGloMid[0]);

				//	compute volume
					bf.Vol = VecTwoNorm(bf.Normal);

				//	compute local integration points
					bf.vWeight = rSCVFQuadRule.weights();
					bf.nip = rSCVFQuadRule.size();
					bf.vLocalIP.resize(bf.nip);
					bf.vGlobalIP.resize(bf.nip);

					ReferenceMapping<scvf_type, dim> map(bf.vLocPos);
					for(size_t ip = 0; ip < rSCVFQuadRule.size(); ++ip)
						map.local_to_global(bf.vLocalIP[ip], rSCVFQuadRule.point(ip));

				//	compute global integration points
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
						rMapping.local_to_global(bf.vGlobalIP[ip], bf.vLocalIP[ip]);

					bf.nsh = rTrialSpace.num_sh();
					bf.vvShape.resize(bf.nip);
					bf.vvLocalGrad.resize(bf.nip);
					bf.vvGlobalGrad.resize(bf.nip);
					bf.vJtInv.resize(bf.nip);
					bf.vDetJ.resize(bf.nip);
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
					{
						bf.vvShape[ip].resize(bf.nsh);
						bf.vvLocalGrad[ip].resize(bf.nsh);
						bf.vvGlobalGrad[ip].resize(bf.nsh);
					}

				//	compute shapes and gradients
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
					{
						rTrialSpace.shapes(&(bf.vvShape[ip][0]), bf.local_ip(ip));
						rTrialSpace.grads(&(bf.vvLocalGrad[ip][0]), bf.local_ip(ip));
					}

					rMapping.jacobian_transposed_inverse(&bf.vJtInv[0], &bf.vLocalIP[0], bf.num_ip());
					rMapping.sqrt_gram_det(&bf.vDetJ[0], &bf.vLocalIP[0], bf.num_ip());

				//	compute global gradient
					for(size_t ip = 0; ip < bf.num_ip(); ++ip)
						for(size_t sh = 0 ; sh < bf.num_sh(); ++sh)
							MatVecMult(bf.vvGlobalGrad[ip][sh],
							           bf.vJtInv[ip], bf.vvLocalGrad[ip][sh]);

				//	increase curr_bf
					++curr_bf;
				}
			}
		}
	}

	}
	UG_CATCH_THROW("DimFVGeometry: update failed.");
}

//////////////////////
// FVGeometry
#ifdef UG_DIM_1
//template class DimFVGeometry<1, 1>;
//template class DimFVGeometry<2, 1>;
//template class DimFVGeometry<3, 1>;
#endif

#ifdef UG_DIM_2
template class FVGeometry<1, Triangle, 2>;
template class FVGeometry<1, Quadrilateral, 2>;
template class FVGeometry<2, Triangle, 2>;
template class FVGeometry<2, Quadrilateral, 2>;
template class FVGeometry<3, Triangle, 2>;
template class FVGeometry<3, Quadrilateral, 2>;

template class DimFVGeometry<2, 2>;
template class DimFVGeometry<3, 2>;
#endif

#ifdef UG_DIM_3
template class FVGeometry<1, Tetrahedron, 3>;
template class FVGeometry<1, Prism, 3>;
template class FVGeometry<1, Hexahedron, 3>;
template class FVGeometry<2, Tetrahedron, 3>;
template class FVGeometry<2, Prism, 3>;
template class FVGeometry<2, Hexahedron, 3>;
template class FVGeometry<3, Tetrahedron, 3>;
template class FVGeometry<3, Prism, 3>;
template class FVGeometry<3, Hexahedron, 3>;

template class DimFVGeometry<3, 3>;
#endif

} // end namespace ug

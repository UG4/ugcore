
/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FV1FT_GEOMETRY_IMPL_
#define __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FV1FT_GEOMETRY_IMPL_

#include "common/util/provider.h"
#include "fv1FT_geom.h"
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

// 	for OCTAHEDRONS: add midpoints of imaginary faces, edges and volumes
// 	resulting from the division into 4 tetrahedra alongside inner edge 3->1
	if (rRefElem.roid() == ROID_OCTAHEDRON)
	{
	// 	diagonal 3->1, diagonal 1->3
		VecScaleAdd(vvMid[1][rRefElem.num(1)], 0.5, vCorner[3], 0.5, vCorner[1]);
		VecScaleAdd(vvMid[1][rRefElem.num(1)+1], 0.5, vCorner[1], 0.5, vCorner[3]);

	// 	subface 1,2,3; subface 1,3,2; subface 1,3,4; subface 1,4,3; face 1,3,5; face 1,5,3; face 1,0,3; face 1,3,0
		vvMid[2][rRefElem.num(2)] = vCorner[1];
		vvMid[2][rRefElem.num(2)] += vCorner[2];
		vvMid[2][rRefElem.num(2)] += vCorner[3];
		vvMid[2][rRefElem.num(2)] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+1] = vCorner[1];
		vvMid[2][rRefElem.num(2)+1] += vCorner[3];
		vvMid[2][rRefElem.num(2)+1] += vCorner[2];
		vvMid[2][rRefElem.num(2)+1] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+2] = vCorner[1];
		vvMid[2][rRefElem.num(2)+2] += vCorner[3];
		vvMid[2][rRefElem.num(2)+2] += vCorner[4];
		vvMid[2][rRefElem.num(2)+2] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+3] = vCorner[1];
		vvMid[2][rRefElem.num(2)+3] += vCorner[4];
		vvMid[2][rRefElem.num(2)+3] += vCorner[3];
		vvMid[2][rRefElem.num(2)+3] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+4] = vCorner[1];
		vvMid[2][rRefElem.num(2)+4] += vCorner[3];
		vvMid[2][rRefElem.num(2)+4] += vCorner[5];
		vvMid[2][rRefElem.num(2)+4] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+5] = vCorner[1];
		vvMid[2][rRefElem.num(2)+5] += vCorner[5];
		vvMid[2][rRefElem.num(2)+5] += vCorner[3];
		vvMid[2][rRefElem.num(2)+5] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+6] = vCorner[1];
		vvMid[2][rRefElem.num(2)+6] += vCorner[0];
		vvMid[2][rRefElem.num(2)+6] += vCorner[3];
		vvMid[2][rRefElem.num(2)+6] *= 1.0/3.0;

		vvMid[2][rRefElem.num(2)+7] = vCorner[1];
		vvMid[2][rRefElem.num(2)+7] += vCorner[3];
		vvMid[2][rRefElem.num(2)+7] += vCorner[0];
		vvMid[2][rRefElem.num(2)+7] *= 1.0/3.0;

	// 	subvolume 1,2,3,5; subvolume 1,3,4,5; subvolume 1,2,3,0; subvolume 1,3,4,0
		vvMid[3][rRefElem.num(3)] = vCorner[1];
		vvMid[3][rRefElem.num(3)] += vCorner[2];
		vvMid[3][rRefElem.num(3)] += vCorner[3];
		vvMid[3][rRefElem.num(3)] += vCorner[5];
		vvMid[3][rRefElem.num(3)] *= 0.25;

		vvMid[3][rRefElem.num(3)+1] = vCorner[1];
		vvMid[3][rRefElem.num(3)+1] += vCorner[3];
		vvMid[3][rRefElem.num(3)+1] += vCorner[4];
		vvMid[3][rRefElem.num(3)+1] += vCorner[5];
		vvMid[3][rRefElem.num(3)+1] *= 0.25;

		vvMid[3][rRefElem.num(3)+2] = vCorner[1];
		vvMid[3][rRefElem.num(3)+2] += vCorner[2];
		vvMid[3][rRefElem.num(3)+2] += vCorner[3];
		vvMid[3][rRefElem.num(3)+2] += vCorner[0];
		vvMid[3][rRefElem.num(3)+2] *= 0.25;

		vvMid[3][rRefElem.num(3)+3] = vCorner[1];
		vvMid[3][rRefElem.num(3)+3] += vCorner[3];
		vvMid[3][rRefElem.num(3)+3] += vCorner[4];
		vvMid[3][rRefElem.num(3)+3] += vCorner[0];
		vvMid[3][rRefElem.num(3)+3] *= 0.25;
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

	if (rRefElem.roid() != ROID_PYRAMID && rRefElem.roid() != ROID_OCTAHEDRON)
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
// 	pyramid here
	else if(rRefElem.roid() == ROID_PYRAMID)
	{
	// 	start at edge midpoint
		vMidID[0] = MidID(1,i/2);

	//  there are 2 scvf per edge
		if(i%2 == 0){
			vMidID[1] = MidID(2, rRefElem.id(1, i/2, 2, 0)); // side 0
			vMidID[2] = MidID(dim, 0); // center of element
		} else {
			vMidID[1] = MidID(dim, 0); // center of element
			vMidID[2] = MidID(2, rRefElem.id(1, i/2, 2, 1)); // side 1
		}
	}
// 	octahedron here (analogue to scvf ordering in pyramids but in top/bottom pairs)
	else if(rRefElem.roid() == ROID_OCTAHEDRON)
	{
		switch (i)
		{
		// 	scvf of edge 4 (top)
			case 0:	vMidID[0] = MidID(1,4);	// edge 4
					vMidID[1] = MidID(2,8);	// subface 0 = 1,2,3
					vMidID[2] = MidID(3,1);	// subvolume 0
					vMidID[3] = MidID(2,4); // face 4
					break;
		// 	scvf of edge 4 (bottom)
			case 1:	vMidID[0] = MidID(1,4);	// edge 4
					vMidID[1] = MidID(2,9);	// subface 1 = 1,3,2
					vMidID[2] = MidID(3,3);	// subvolume 2
					vMidID[3] = MidID(2,0); // face 0
					break;


		// 	scvf of edge 5 (top)
			case 2:	vMidID[0] = MidID(1,5);	// edge 5
					vMidID[1] = MidID(2,8);	// subface 0 = 1,2,3
					vMidID[2] = MidID(3,1);	// subvolume 0
					vMidID[3] = MidID(2,5); // face 5
					break;
		// 	scvf of edge 5 (bottom)
			case 3:	vMidID[0] = MidID(1,5);	// edge 5
					vMidID[1] = MidID(2,9);	// subface 1 = 1,3,2
					vMidID[2] = MidID(3,3);	// subvolume 2
					vMidID[3] = MidID(2,1); // face 1
					break;


		// 	scvf of diagonal 3->1 (top) in subvolume 0
			case 4:	vMidID[0] = MidID(1,12);// diagonal 3->1
					vMidID[1] = MidID(2,8);	// subface 0 = 1,2,3
					vMidID[2] = MidID(3,1);	// subvolume 0
					vMidID[3] = MidID(2,13);// face 1,5,3
					break;
		// 	scvf of diagonal 1->3 (bottom) in subvolume 2
			case 5:	vMidID[0] = MidID(1,13);// diagonal 1->3
					vMidID[1] = MidID(2,9);	// subface 1 = 1,3,2
					vMidID[2] = MidID(3,3);	// subvolume 2
					vMidID[3] = MidID(2,14);// face 1,0,3
					break;


		//	scvf of edge 8 in subvolume 0
			case 6:	vMidID[0] = MidID(1,8);	// edge 8
					vMidID[1] = MidID(2,4); // face 4
					vMidID[2] = MidID(3,1);	// subvolume 0
					vMidID[3] = MidID(2,13);// face 1,5,3
					break;
		//	scvf of edge 0 in subvolume 2
			case 7:	vMidID[0] = MidID(1,0);	// edge 0
					vMidID[1] = MidID(2,15);// face 1,3,0
					vMidID[2] = MidID(3,3);	// subvolume 2
					vMidID[3] = MidID(2,0); // face 0
					break;


		// 	scvf of edge 9
			case 8:	vMidID[0] = MidID(1,9);	// edge 9
					vMidID[1] = MidID(2,5);	// face 5
					vMidID[2] = MidID(3,1);	// subvolume 0
					vMidID[3] = MidID(2,4); // face 4
					break;
		// 	scvf of edge 1
			case 9:	vMidID[0] = MidID(1,1);	// edge 1
					vMidID[1] = MidID(2,0); // face 0
					vMidID[2] = MidID(3,3);	// subvolume 2
					vMidID[3] = MidID(2,1);	// face 1
					break;


		// 	scvf of edge 10 in subvolume 0
			case 10:vMidID[0] = MidID(1,10);// edge 10
					vMidID[1] = MidID(2,12);// face 1,3,5
					vMidID[2] = MidID(3,1);	// subvolume 0
					vMidID[3] = MidID(2,5);	// face 5
					break;
		// 	scvf of edge 2 in subvolume 2
			case 11:vMidID[0] = MidID(1,2);	// edge 2
					vMidID[1] = MidID(2,1); // face 1
					vMidID[2] = MidID(3,3);	// subvolume 2
					vMidID[3] = MidID(2,14);// face 1,0,3
					break;


		// 	scvf of diagonal 1->3 (top) in subvolume 1
			case 12:vMidID[0] = MidID(1,13);// diagonal 1->3
					vMidID[1] = MidID(2,10);// subface 2 = 1,3,4
					vMidID[2] = MidID(3,2);	// subvolume 1
					vMidID[3] = MidID(2,12);// face 1,3,5
					break;
		// 	scvf of diagonal 3->1 (bottom) in subvolume 3
			case 13:vMidID[0] = MidID(1,12);// diagonal 3->1
					vMidID[1] = MidID(2,11);// subface 3 = 1,4,3
					vMidID[2] = MidID(3,4);	// subvolume 3
					vMidID[3] = MidID(2,15);// face 1,3,0
					break;


		// 	scvf of edge 6 (top)
			case 14:vMidID[0] = MidID(1,6);	// edge 6
					vMidID[1] = MidID(2,10);// subface 2 = 1,3,4
					vMidID[2] = MidID(3,2);	// subvolume 1
					vMidID[3] = MidID(2,6); // face 6
					break;
		// 	scvf of edge 6 (bottom)
			case 15:vMidID[0] = MidID(1,6);	// edge 6
					vMidID[1] = MidID(2,11);// subface 3 = 1,4,3
					vMidID[2] = MidID(3,4);	// subvolume 3
					vMidID[3] = MidID(2,2);	// face 2
					break;


		// 	scvf of edge 7 (top)
			case 16:vMidID[0] = MidID(1,7);	// edge 7
					vMidID[1] = MidID(2,10);// subface 2 = 1,3,4
					vMidID[2] = MidID(3,2);	// subvolume 1
					vMidID[3] = MidID(2,7); // face 7
					break;
		// 	scvf of edge 7 (bottom)
			case 17:vMidID[0] = MidID(1,7);	// edge 7
					vMidID[1] = MidID(2,11);// subface 3 = 1,4,3
					vMidID[2] = MidID(3,4);	// subvolume 3
					vMidID[3] = MidID(2,3);	// face 3
					break;


		// 	scvf of edge 8 in subvolume 1
			case 18:vMidID[0] = MidID(1,8);	// edge 8
					vMidID[1] = MidID(2,13);// face 1,5,3
					vMidID[2] = MidID(3,2);	// subvolume 1
					vMidID[3] = MidID(2,7); // face 7
					break;
		// 	scvf of edge 0 in subvolume 3
			case 19:vMidID[0] = MidID(1,0);	// edge 0
					vMidID[1] = MidID(2,3);	// face 3
					vMidID[2] = MidID(3,4);	// subvolume 3
					vMidID[3] = MidID(2,15);// face 1,3,0
					break;


		// 	scvf of edge 10 in subvolume 1
			case 20:vMidID[0] = MidID(1,10);// edge 10
					vMidID[1] = MidID(2,6); // face 6
					vMidID[2] = MidID(3,2);	// subvolume 1
					vMidID[3] = MidID(2,12);// face 1,3,5
					break;
		// 	scvf of edge 2 in subvolume 3
			case 21:vMidID[0] = MidID(1,2);	// edge 2
					vMidID[1] = MidID(2,14);// face 1,0,3
					vMidID[2] = MidID(3,4);	// subvolume 3
					vMidID[3] = MidID(2,2); // face 2
					break;


		// 	scvf of edge 11 in subvolume 1
			case 22:vMidID[0] = MidID(1,11);// edge 11
					vMidID[1] = MidID(2,7); // face 7
					vMidID[2] = MidID(3,2);	// subvolume 1
					vMidID[3] = MidID(2,6); // face 6
					break;
		// 	scvf of edge 3 in subvolume 3
			case 23:vMidID[0] = MidID(1,3);	// edge 3
					vMidID[1] = MidID(2,2); // face 2
					vMidID[2] = MidID(3,4);	// subvolume 3
					vMidID[3] = MidID(2,3); // face 3
					break;


			default:UG_THROW("Octahedron only has 24 SCVFs (no. 0-23), but requested no. " << i << ".");
					break;
		}
	}
}

/**
 * \param[in]	i		indicates that scv corresponds to i'th corner of ref elem
 */
template <typename TRefElem>
static void ComputeSCVMidID(const TRefElem& rRefElem,
                            MidID vMidID[], int i)
{
	static const int dim = TRefElem::dim;

	if (rRefElem.roid() != ROID_PYRAMID && rRefElem.roid() != ROID_OCTAHEDRON)
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
		else {UG_THROW("Dimension higher than 3 not implemented.");}
	}
// 	pyramid here
	else if (rRefElem.roid() == ROID_PYRAMID)
	{
	// 	start at edge midpoint
		vMidID[3] = MidID(1,i/4);

	//  there are 2 scvf per edge
		if(i%4 == 0 || i%4 == 1){
			vMidID[1] = MidID(2, rRefElem.id(1, i/4, 2, 0)); // side 0
			vMidID[2] = MidID(dim, 0); // center of element
		} else {
			vMidID[1] = MidID(dim, 0); // center of element
			vMidID[2] = MidID(2, rRefElem.id(1, i/4, 2, 1)); // side 1
		}

	//	connect to from / to corners of edge
		if(i%2 == 0){
			vMidID[0] = MidID(0, rRefElem.id(1, i/4, 0, 0)); // from
		} else {
			vMidID[0] = MidID(0, rRefElem.id(1, i/4, 0, 1)); // to
		}
	}
	// octahedron here (analogue to scv ordering in pyramids but in top/bottom pairs)
	else if(rRefElem.roid() == ROID_OCTAHEDRON)
	{
		switch (i)
		{
		// 	scv of corner 1 in subvolume 0 (top)
			case 0:	vMidID[0] = MidID(0,1);	// corner 1
					vMidID[1] = MidID(1,4);	// edge 4
					vMidID[2] = MidID(2,8);	// subface 0 = 1,2,3
					vMidID[3] = MidID(1,12);// edge 3->1
					vMidID[4] = MidID(1,8);	// edge 8
					vMidID[5] = MidID(2,4); // face 4
					vMidID[6] = MidID(3,1);	// subvolume 0
					vMidID[7] = MidID(2,13);// face 1,5,3
					break;
		// 	scv of corner 1 in subvolume 2 (bottom)
			case 1:	vMidID[0] = MidID(0,1);	// corner 1
					vMidID[1] = MidID(1,13);// edge 1->3
					vMidID[2] = MidID(2,9);	// subface 1 = 1,3,2
					vMidID[3] = MidID(1,4);	// edge 4
					vMidID[4] = MidID(1,0);	// edge 0
					vMidID[5] = MidID(2,15);// face 1,3,0
					vMidID[6] = MidID(3,3);	// subvolume 2
					vMidID[7] = MidID(2,0); // face 0
					break;


		// 	scv of corner 2 in subvolume 0 (top)
			case 2:	vMidID[0] = MidID(0,2);	// corner 2
					vMidID[1] = MidID(1,5);	// edge 5
					vMidID[2] = MidID(2,8);	// subface 0 = 1,2,3
					vMidID[3] = MidID(1,4); // edge 4
					vMidID[4] = MidID(1,9);	// edge 9
					vMidID[5] = MidID(2,5); // face 5
					vMidID[6] = MidID(3,1);	// subvolume 0
					vMidID[7] = MidID(2,4); // face 4
					break;
		// 	scv of corner 2 in subvolume 2 (bottom)
			case 3:	vMidID[0] = MidID(0,2);	// corner 2
					vMidID[1] = MidID(1,4); // edge 4
					vMidID[2] = MidID(2,9);	// subface 1 = 1,3,2
					vMidID[3] = MidID(1,5);	// edge 5
					vMidID[4] = MidID(1,1);	// edge 1
					vMidID[5] = MidID(2,0); // face 0
					vMidID[6] = MidID(3,3);	// subvolume 2
					vMidID[7] = MidID(2,1); // face 1
					break;


		// 	scv of corner 3 in subvolume 0 (top)
			case 4:	vMidID[0] = MidID(0,3);	// corner 3
					vMidID[1] = MidID(1,12);// edge 3->1
					vMidID[2] = MidID(2,8);	// subface 0 = 1,2,3
					vMidID[3] = MidID(1,5);	// edge 5
					vMidID[4] = MidID(1,10);// edge 10
					vMidID[5] = MidID(2,13);// face 1,5,3
					vMidID[6] = MidID(3,1);	// subvolume 0
					vMidID[7] = MidID(2,5);	// face 5
					break;
		// 	scv of corner 3 in subvolume 2 (bottom)
			case 5:	vMidID[0] = MidID(0,3);	// corner 3
					vMidID[1] = MidID(1,5);	// edge 5
					vMidID[2] = MidID(2,9);	// subface 0 = 1,3,2
					vMidID[3] = MidID(1,13);// edge 1->3
					vMidID[4] = MidID(1,2); // edge 2
					vMidID[5] = MidID(2,1);	// face 1
					vMidID[6] = MidID(3,3);	// subvolume 2
					vMidID[7] = MidID(2,15);// face 1,3,0
					break;


		// 	scv of corner 5 in subvolume 0 (top)
			case 6:	vMidID[0] = MidID(0,5);	// corner 5
					vMidID[1] = MidID(1,9); // edge 9
					vMidID[2] = MidID(2,4);	// face 4
					vMidID[3] = MidID(1,8);	// edge 8
					vMidID[4] = MidID(1,10);// edge 10
					vMidID[5] = MidID(2,5); // face 5
					vMidID[6] = MidID(3,1);	// subvolume 0
					vMidID[7] = MidID(2,13);// subface 1,5,3
					break;
		// 	scv of corner 0 in subvolume 2 (bottom)
			case 7:	vMidID[0] = MidID(0,0);	// corner 0
					vMidID[1] = MidID(1,0);	// edge 0
					vMidID[2] = MidID(2,0);	// face 0
					vMidID[3] = MidID(1,1); // edge 1
					vMidID[4] = MidID(1,2); // edge 2
					vMidID[5] = MidID(2,14);// subface 1,0,3
					vMidID[6] = MidID(3,3);	// subvolume 2
					vMidID[7] = MidID(2,1); // face 1
					break;


		// 	scv of corner 1 in subvolume 1 (top)
			case 8:	vMidID[0] = MidID(0,1);	// corner 1
					vMidID[1] = MidID(1,13);// edge 1->3
					vMidID[2] = MidID(2,10);// subface 2 = 1,3,4
					vMidID[3] = MidID(1,7); // edge 7
					vMidID[4] = MidID(1,8);	// edge 8
					vMidID[5] = MidID(2,12);// face 1,3,5
					vMidID[6] = MidID(3,2);	// subvolume 1
					vMidID[7] = MidID(2,7); // face 7
					break;
		//	scv of corner 1 in subvolume 3 (bottom)
			case 9:	vMidID[0] = MidID(0,1);	// corner 1
					vMidID[1] = MidID(1,7); // edge 7
					vMidID[2] = MidID(2,11);// subface 3 = 1,4,3
					vMidID[3] = MidID(1,12);// edge 3->1
					vMidID[4] = MidID(1,0);	// edge 0
					vMidID[5] = MidID(2,3); // face 3
					vMidID[6] = MidID(3,4);	// subvolume 3
					vMidID[7] = MidID(2,14);// face 1,0,3
					break;


		// 	scv of corner 3 in subvolume 1 (top)
			case 10:vMidID[0] = MidID(0,3);	// corner 3
					vMidID[1] = MidID(1,6); // edge 6
					vMidID[2] = MidID(2,10);// subface 2 = 1,3,4
					vMidID[3] = MidID(1,13);// edge 1->3
					vMidID[4] = MidID(1,10);// edge 10
					vMidID[5] = MidID(2,6); // face 6
					vMidID[6] = MidID(3,2);	// subvolume 1
					vMidID[7] = MidID(2,12);// face 1,3,5
					break;
		// 	scv of corner 3 in subvolume 3 (bottom)
			case 11:vMidID[0] = MidID(0,3);	// corner 3
					vMidID[1] = MidID(1,12);// edge 3->1
					vMidID[2] = MidID(2,11);// subface 3 = 1,4,3
					vMidID[3] = MidID(1,6); // edge 6
					vMidID[4] = MidID(1,2);	// edge 2
					vMidID[5] = MidID(2,14);// face 1,0,3
					vMidID[6] = MidID(3,4);	// subvolume 3
					vMidID[7] = MidID(2,2); // face 2
					break;


		// 	scv of corner 4 in subvolume 1 (top)
			case 12:vMidID[0] = MidID(0,4);	// corner 4
					vMidID[1] = MidID(1,7); // edge 7
					vMidID[2] = MidID(2,10);// subface 2 = 1,3,4
					vMidID[3] = MidID(1,6); // edge 6
					vMidID[4] = MidID(1,11);// edge 11
					vMidID[5] = MidID(2,7); // face 7
					vMidID[6] = MidID(3,2);	// subvolume 1
					vMidID[7] = MidID(2,6); // face 6
					break;
		// 	scv of corner 4 in subvolume 3 (bottom)
			case 13:vMidID[0] = MidID(0,4);	// corner 4
					vMidID[1] = MidID(1,6); // edge 6
					vMidID[2] = MidID(2,11);// subface 3 = 1,4,3
					vMidID[3] = MidID(1,7); // edge 7
					vMidID[4] = MidID(1,3);	// edge 3
					vMidID[5] = MidID(2,2); // face 2
					vMidID[6] = MidID(3,4);	// subvolume 3
					vMidID[7] = MidID(2,3); // face 3
					break;


		// 	scv of corner 5 in subvolume 1 (top)
			case 14:vMidID[0] = MidID(0,5);	// corner 5
					vMidID[1] = MidID(1,10);// edge 10
					vMidID[2] = MidID(2,12);// subface 1,3,5
					vMidID[3] = MidID(1,8);	// edge 8
					vMidID[4] = MidID(1,11);// edge 11
					vMidID[5] = MidID(2,6);	// face 6
					vMidID[6] = MidID(3,2);	// subvolume 1
					vMidID[7] = MidID(2,7); // face 7
					break;
		// 	scv of corner 0 in subvolume 3 (bottom)
			case 15:vMidID[0] = MidID(0,0);	// corner 0
					vMidID[1] = MidID(1,0);	// edge 0
					vMidID[2] = MidID(2,15);// subface 1,3,0
					vMidID[3] = MidID(1,2); // edge 2
					vMidID[4] = MidID(1,3); // edge 3
					vMidID[5] = MidID(2,3); // face 3
					vMidID[6] = MidID(3,4);	// subvolume 3
					vMidID[7] = MidID(2,2); // face 2
					break;
			default:UG_THROW("Octahedron only has 16 SCVs (no. 0-15), but requested no. " << i << ".");
					break;
		}
	}
}

/**
 * \param[in]	i		indicates that bf corresponds to i'th corner of ref elem
 */
template <typename TRefElem>
static void ComputeBFMidID(const TRefElem& rRefElem, int side,
                            MidID vMidID[], int co)
{
	static const int dim = TRefElem::dim;

	//	number of corners of side
	const int coOfSide = rRefElem.num(dim-1, side, 0);

	// 	set mid ids
	if (dim == 1)
	{
		vMidID[0] = MidID(0, rRefElem.id(0, side, 0, co));
	}
	else if(dim == 2)
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
// FV1FT Geometry for Reference Element Type
////////////////////////////////////////////////////////////////////////////////

template <typename TElem, int TWorldDim>
FV1FTGeometry<TElem, TWorldDim>::
FV1FTGeometry()
	: m_pElem(NULL), m_rRefElem(Provider<ref_elem_type>::get()),
	  m_rTrialSpace(Provider<local_shape_fct_set_type>::get())
{
	update_local_data();
}

template <typename TElem, int TWorldDim>
void FV1FTGeometry<TElem, TWorldDim>::
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
		if (m_rRefElem.REFERENCE_OBJECT_ID != ROID_PYRAMID && m_rRefElem.REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
		{
			m_vSCVF[i].From = m_rRefElem.id(1, i, 0, 0);
			m_vSCVF[i].To = m_rRefElem.id(1, i, 0, 1);
		}
		// special case pyramid
		else if(m_rRefElem.REFERENCE_OBJECT_ID == ROID_PYRAMID)
		{
		// 	map according to order defined in ComputeSCVFMidID
			m_vSCVF[i].From = m_rRefElem.id(1, i/2, 0, 0);
			m_vSCVF[i].To = m_rRefElem.id(1, i/2, 0, 1);
		}
		//	special case octahedron (scvf not mappable by edges)
		else if(m_rRefElem.REFERENCE_OBJECT_ID == ROID_OCTAHEDRON)
		{
		// 	map according to order defined in ComputeSCVFMidID
			switch(i)
			{
				case 0:	m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 2;
						break;
				case 1:	m_vSCVF[i].From = 2;
						m_vSCVF[i].To	= 1;
						break;


				case 2:	m_vSCVF[i].From = 2;
						m_vSCVF[i].To	= 3;
						break;
				case 3:	m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 2;
						break;


				case 4:	m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 1;
						break;
				case 5:	m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 3;
						break;


				case 6:	m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 5;
						break;
				case 7:	m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 0;
						break;


				case 8:	m_vSCVF[i].From = 2;
						m_vSCVF[i].To	= 5;
						break;
				case 9:	m_vSCVF[i].From = 2;
						m_vSCVF[i].To	= 0;
						break;


				case 10:m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 5;
						break;
				case 11:m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 0;
						break;


				case 12:m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 3;
						break;
				case 13:m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 1;
						break;


				case 14:m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 4;
						break;
				case 15:m_vSCVF[i].From = 4;
						m_vSCVF[i].To	= 3;
						break;


				case 16:m_vSCVF[i].From = 4;
						m_vSCVF[i].To	= 1;
						break;
				case 17:m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 4;
						break;


				case 18:m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 5;
						break;
				case 19:m_vSCVF[i].From = 1;
						m_vSCVF[i].To	= 0;
						break;


				case 20:m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 5;
						break;
				case 21:m_vSCVF[i].From = 3;
						m_vSCVF[i].To	= 0;
						break;


				case 22:m_vSCVF[i].From = 4;
						m_vSCVF[i].To	= 5;
						break;
				case 23:m_vSCVF[i].From = 4;
						m_vSCVF[i].To	= 0;
						break;
			}
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
		if (m_rRefElem.REFERENCE_OBJECT_ID != ROID_PYRAMID && m_rRefElem.REFERENCE_OBJECT_ID != ROID_OCTAHEDRON)
		{
			m_vSCV[i].nodeId = i;
		}
		// special case pyramid
		else if(m_rRefElem.REFERENCE_OBJECT_ID == ROID_PYRAMID)
		{
		// 	map according to order defined in ComputeSCVMidID
			if(i%2 == 0){
				m_vSCV[i].nodeId  = m_rRefElem.id(1, i/4, 0, 0); // from
			} else {
				m_vSCV[i].nodeId  = m_rRefElem.id(1, i/4, 0, 1); // to
			}
		}
		// special case octahedron (scvf not mappable by edges)
		else if(m_rRefElem.REFERENCE_OBJECT_ID == ROID_OCTAHEDRON)
		{
		// 	map according to order defined in ComputeSCVMidID
			switch(i)
			{
				case 0: m_vSCV[i].nodeId = 1;
						break;
				case 1: m_vSCV[i].nodeId = 1;
						break;


				case 2: m_vSCV[i].nodeId = 2;
						break;
				case 3: m_vSCV[i].nodeId = 2;
						break;


				case 4: m_vSCV[i].nodeId = 3;
						break;
				case 5: m_vSCV[i].nodeId = 3;
						break;


				case 6: m_vSCV[i].nodeId = 5;
						break;
				case 7: m_vSCV[i].nodeId = 0;
						break;


				case 8: m_vSCV[i].nodeId = 1;
						break;
				case 9: m_vSCV[i].nodeId = 1;
						break;


				case 10:m_vSCV[i].nodeId = 3;
						break;
				case 11:m_vSCV[i].nodeId = 3;
						break;


				case 12:m_vSCV[i].nodeId = 4;
						break;
				case 13:m_vSCV[i].nodeId = 4;
						break;


				case 14:m_vSCV[i].nodeId = 5;
						break;
				case 15:m_vSCV[i].nodeId = 0;
						break;
			}
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

	if(ref_elem_type::REFERENCE_OBJECT_ID == ROID_PYRAMID || ref_elem_type::REFERENCE_OBJECT_ID == ROID_OCTAHEDRON)
		for(size_t i = 0; i < num_scv(); ++i)
			m_vLocSCV_IP[i] = scv(i).local_ip();
}

/// update data for given element
template <typename TElem, int TWorldDim>
void FV1FTGeometry<TElem, TWorldDim>::
update(GridObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

// 	if already update for this element, do nothing
	if(m_pElem == pElem) return; else m_pElem = pElem;

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

		/*
		 *	Only for debug purposes testing octahedral FV1 discretization
		 *
			MathVector<dim> baryCenter(0.0);
			for(size_t j = 0; j < 8; ++j)
			{
				baryCenter += m_vSCV[i].vLocPos[j];
			}

			baryCenter *= 1.0/8.0;
		*
		*/

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

	if(ref_elem_type::REFERENCE_OBJECT_ID == ROID_PYRAMID || ref_elem_type::REFERENCE_OBJECT_ID == ROID_OCTAHEDRON)
		for(size_t i = 0; i < num_scv(); ++i)
			m_vGlobSCV_IP[i] = scv(i).global_ip();

//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}

template <typename TElem, int TWorldDim>
void FV1FTGeometry<TElem, TWorldDim>::
update_boundary_faces(GridObject* elem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
	UG_ASSERT(dynamic_cast<TElem*>(elem) != NULL, "Wrong element type.");
	TElem* pElem = static_cast<TElem*>(elem);

//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<Vertex*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<Edge*> vEdges;
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

		//	number of corners of side
			const int coOfSide = m_rRefElem.num(dim-1, side, 0);

		//	resize vector
			vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				BF& bf = vBF[curr_bf];

			//	set node id == scv this bf belongs to
				bf.nodeId = m_rRefElem.id(dim-1, side, 0, co);

			//	Compute MidID for BF
				ComputeBFMidID(m_rRefElem, side, bf.vMidID, co);

			// 	copy corners of bf
				CopyCornerByMidID<dim, maxMid>(bf.vLocPos, bf.vMidID, m_vvLocMid, BF::numCo);
				CopyCornerByMidID<worldDim, maxMid>(bf.vGloPos, bf.vMidID, m_vvGloMid, BF::numCo);

			// 	integration point
				AveragePositions(bf.localIP, bf.vLocPos, BF::numCo);
				AveragePositions(bf.globalIP, bf.vGloPos, BF::numCo);

			// 	normal on scvf
				traits::NormalOnBF(bf.Normal, bf.vGloPos, m_vvGloMid[0]);

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
template <int TDim, int TWorldDim, class TInterfaceHandler>
void DimFV1FTGeometry<TDim, TWorldDim,  TInterfaceHandler>::
update_local(ReferenceObjectID roid)
{
 	m_roid = roid;
 
//	get reference element
	try
	{
 		const DimReferenceElement<dim>& rRefElem
			= ReferenceElementProvider::get<dim>(m_roid);
	
 	// 	set corners of element as local centers of nodes
		for(size_t i = 0; i < rRefElem.num(0); ++i)
			m_vvLocMid[0][i] = rRefElem.corner(i);
	
	//	compute local midpoints
		ComputeMidPoints<dim, DimReferenceElement<dim>, maxMid>
											(rRefElem, m_vvLocMid[0], m_vvLocMid);
	
	//	set number of scvf / scv of this roid
		if(m_roid != ROID_PYRAMID && m_roid != ROID_OCTAHEDRON)
		{
			m_numSCV  = rRefElem.num(0);
			m_numSCVF = rRefElem.num(1);
		}
		else if(dim == 3 && m_roid == ROID_PYRAMID)
		{
			UG_WARNING("Pyramid Finite Volume Geometry for 1st order currently  "
					"implemented in DimFV1Geom EXPERIMENTATLLY. Please contact "
					"Martin Stepniewski or Andreas Vogel if you see this message.")

			m_numSCV  = 4*rRefElem.num(1);
			m_numSCVF = 2*rRefElem.num(1);
		}
		else if(dim == 3 && m_roid == ROID_OCTAHEDRON)
		{
		// 	Case octahedron
			m_numSCV  = 16;
			m_numSCVF = 24;
		}
	
	// 	set up local informations for SubControlVolumeFaces (scvf)
	// 	each scvf is associated to one edge of the element
		for(size_t i = 0; i < num_scvf(); ++i)
		{
		//	this scvf separates the given nodes
			if (m_roid != ROID_PYRAMID && m_roid != ROID_OCTAHEDRON)
			{
				m_vSCVF[i].From = rRefElem.id(1, i, 0, 0);
				m_vSCVF[i].To = rRefElem.id(1, i, 0, 1);
			}
			// special case pyramid (scvf not mappable by edges)
			else if (dim == 3 && m_roid == ROID_PYRAMID)
			{
			// 	map according to order defined in ComputeSCVFMidID
				m_vSCVF[i].From = rRefElem.id(1, i/2, 0, 0);
				m_vSCVF[i].To = rRefElem.id(1, i/2, 0, 1);
			}
			//	special case octahedron (scvf not mappable by edges)
			else if(dim == 3 && m_roid == ROID_OCTAHEDRON)
			{
			// 	map according to order defined in ComputeSCVFMidID
				switch(i)
				{
					case 0:	m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 2;
							break;
					case 1:	m_vSCVF[i].From = 2;
							m_vSCVF[i].To	= 1;
							break;


					case 2:	m_vSCVF[i].From = 2;
							m_vSCVF[i].To	= 3;
							break;
					case 3:	m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 2;
							break;


					case 4:	m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 1;
							break;
					case 5:	m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 3;
							break;


					case 6:	m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 5;
							break;
					case 7:	m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 0;
							break;


					case 8:	m_vSCVF[i].From = 2;
							m_vSCVF[i].To	= 5;
							break;
					case 9:	m_vSCVF[i].From = 2;
							m_vSCVF[i].To	= 0;
							break;


					case 10:m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 5;
							break;
					case 11:m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 0;
							break;


					case 12:m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 3;
							break;
					case 13:m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 1;
							break;


					case 14:m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 4;
							break;
					case 15:m_vSCVF[i].From = 4;
							m_vSCVF[i].To	= 3;
							break;


					case 16:m_vSCVF[i].From = 4;
							m_vSCVF[i].To	= 1;
							break;
					case 17:m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 4;
							break;


					case 18:m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 5;
							break;
					case 19:m_vSCVF[i].From = 1;
							m_vSCVF[i].To	= 0;
							break;


					case 20:m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 5;
							break;
					case 21:m_vSCVF[i].From = 3;
							m_vSCVF[i].To	= 0;
							break;


					case 22:m_vSCVF[i].From = 4;
							m_vSCVF[i].To	= 5;
							break;
					case 23:m_vSCVF[i].From = 4;
							m_vSCVF[i].To	= 0;
							break;
				}
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
			if (m_roid != ROID_PYRAMID && m_roid != ROID_OCTAHEDRON)
			{
				m_vSCV[i].nodeId = i;
			}
			// special case pyramid (scv not mappable by corners)
			else if(dim == 3 && m_roid == ROID_PYRAMID)
			{
			// 	map according to order defined in ComputeSCVMidID
				if(i%2 == 0){
					m_vSCV[i].nodeId  = rRefElem.id(1, i/4, 0, 0); // from
				} else {
					m_vSCV[i].nodeId  = rRefElem.id(1, i/4, 0, 1); // to
				}
			}
			// special case octahedron (scvf not mappable by edges)
			else if(dim == 3 && m_roid == ROID_OCTAHEDRON)
			{
			// 	map according to order defined in ComputeSCVMidID
				switch(i)
				{
					case 0: m_vSCV[i].nodeId = 1;
							break;
					case 1: m_vSCV[i].nodeId = 1;
							break;


					case 2: m_vSCV[i].nodeId = 2;
							break;
					case 3: m_vSCV[i].nodeId = 2;
							break;


					case 4: m_vSCV[i].nodeId = 3;
							break;
					case 5: m_vSCV[i].nodeId = 3;
							break;


					case 6: m_vSCV[i].nodeId = 5;
							break;
					case 7: m_vSCV[i].nodeId = 0;
							break;


					case 8: m_vSCV[i].nodeId = 1;
							break;
					case 9: m_vSCV[i].nodeId = 1;
							break;


					case 10:m_vSCV[i].nodeId = 3;
							break;
					case 11:m_vSCV[i].nodeId = 3;
							break;


					case 12:m_vSCV[i].nodeId = 4;
							break;
					case 13:m_vSCV[i].nodeId = 4;
							break;


					case 14:m_vSCV[i].nodeId = 5;
							break;
					case 15:m_vSCV[i].nodeId = 0;
							break;
				}
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
			LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, dim, 1));
	
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
	UG_CATCH_THROW("DimFV1FTGeometry: update failed.");

// 	copy ip positions in a list for Sub Control Volumes Faces (SCVF)
	for(size_t i = 0; i < num_scvf(); ++i)
		m_vLocSCVF_IP[i] = scvf(i).local_ip();

	if(roid == ROID_PYRAMID || roid == ROID_OCTAHEDRON)
		for(size_t i = 0; i < num_scv(); ++i)
			m_vLocSCV_IP[i] = scv(i).local_ip();
}


/// update data for given element
template <int TDim, int TWorldDim, class TInterfaceHandler>
void DimFV1FTGeometry<TDim, TWorldDim,  TInterfaceHandler>::
update(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
  
///////////////////////////////////////////////////////////////////////////////////
// (1) - (4= pre-processing for re-computing the cooridnates of element
// corners 'vCornerCoords' in case of a cut element
    
    set_element_modus(false);
    
// If already update for this element, do nothing
   if(m_pElem == pElem) return; else m_pElem = pElem;
   
// (1) collect potentially new 'vCornerCoords' and set the elemModus
//       (INSIDE_DOM/OUTSIDE_DOM/CUT_BY_INTERFACE)
    bool do_update_local = false;
    try{
        do_update_local = m_spInterfaceHandler->update_elem(pElem, vCornerCoords);
    }
    UG_CATCH_THROW("DimFV1FTGeometry: access to InterfaceHander::update_elem() failed!.");

    
// (2) If OUTSIDE_DOM: set m_roid to usual value and perform local update
   if( m_spInterfaceHandler->elementModus() == OUTSIDE_DOM)
   {
       if ( dim == 2 )      m_roid = ROID_TRIANGLE;
       if ( dim == 3 )      m_roid = ROID_TETRAHEDRON;

       update_local(m_roid);
   }
   
// (3) get computed roid and perform the local update based on that
   if (do_update_local || m_roid != m_spInterfaceHandler->roid() )
   {
       if( m_spInterfaceHandler->elementModus() != OUTSIDE_DOM )
       {
           m_roid = m_spInterfaceHandler->roid();
           update_local(m_roid);
       }
   }
    
// (4) set needed flag
   if(m_spInterfaceHandler->elementModus() == CUT_BY_INTERFACE)
       set_element_modus(true);
   
    
//////////////////////////////////////////////////////////////////////////////////////
// now: perform usual computations on potentially updated cut element 'vCornerCoords'

//	get reference element
	try{
	const DimReferenceElement<dim>& rRefElem
		= ReferenceElementProvider::get<dim>(m_roid);

// 	remember global position of nodes
	for(size_t i = 0; i < rRefElem.num(0); ++i)
        m_vvGloMid[0][i] = m_spInterfaceHandler->corner(i);
 
        
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
    rMapping.update(m_spInterfaceHandler->corners());
    //rMapping.update(vCornerCoords);


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

	if(m_roid == ROID_PYRAMID || m_roid == ROID_OCTAHEDRON)
		for(size_t i = 0; i < num_scv(); ++i)
			m_vGlobSCV_IP[i] = scv(i).global_ip();

	}
	UG_CATCH_THROW("DimFV1FTGeometry: update failed.");


/////////////////////////////////////////////////////////////////////////////
// compute boundary faces of inner boundary: m_vBF-data
//      --> needed for assembling boundary-conditions on immersed interface
    if(m_spInterfaceHandler->elementModus() == CUT_BY_INTERFACE)
    {
        update_inner_boundary_faces();
        m_spInterfaceHandler->update_inner_boundary(m_spInterfaceHandler->corners());

    }
    else if(m_spInterfaceHandler->elementModus() == CUT_BY_2_INTERFACE)
    {
        if ( m_spInterfaceHandler->StdFV_assembling() )
        {
            update_inner_boundary_faces();
            m_spInterfaceHandler->update_inner_boundary_for2();
        }
        else
        {
            update_inner_boundary_faces_for2();
            m_spInterfaceHandler->update_inner_boundary_for2();
        }
    }

     
    
//	if no boundary subsets required, return
	if(num_boundary_subsets() == 0 || ish == NULL) return;
	else update_boundary_faces(pElem, vCornerCoords, ish);
}
    
    
    
// compare implementation of 'DimFV1Geometry::update_boundary_faces()'
template <int TDim, int TWorldDim, class TInterfaceHandler>
void DimFV1FTGeometry<TDim, TWorldDim,  TInterfaceHandler>::
update_inner_boundary_faces()
{
    /////////////////////////////////////////////////////////////////////////////
    // get general data
    /////////////////////////////////////////////////////////////////////////////
    try{
        const DimReferenceElement<dim>& rRefElem
        = ReferenceElementProvider::get<dim>(m_roid);
        
        DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(m_roid);
        rMapping.update(m_spInterfaceHandler->corners());
            
        const LocalShapeFunctionSet<dim>& TrialSpace =
        LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, dim, 1));
        
        /////////////////////////////////////////////////////////////////////////////
        // collect boudary faces
        /////////////////////////////////////////////////////////////////////////////
            
        // get number of sides of element
        size_t numSides = 0;
        numSides = rRefElem.num(dim-1);
        
        //	current number of bf
        size_t curr_bf = 0;
        
        //	std::vector<BF>& vBF = m_spInterfaceHandler->m_vBF;
        std::vector<BF>& vBF = m_spInterfaceHandler->get_boundary_faces();
 
        vBF.clear();
        
        //	loop sides of element
        for(size_t side = 0; side < numSides; ++side)
        {
            // side is no boundary face => continue
            if( m_spInterfaceHandler->elementModus() == CUT_BY_2_INTERFACE)
            {
                if ( !m_spInterfaceHandler->is_boundary_face_for2(side) )
                    continue;
            }
            else
            {
                if ( !m_spInterfaceHandler->is_boundary_face(side) )
                    continue;
            }
            //	number of corners of side (special case bottom side pyramid)
            const int coOfSide = (m_roid != ROID_PYRAMID || side != 0)
            ? rRefElem.num(dim-1, side, 0) : rRefElem.num(dim-1, side, 0) + 2;
            
            //	resize vector
            vBF.resize(curr_bf + coOfSide);
                
            //	loop corners
            for(int co = 0; co < coOfSide; ++co)
            {
                //	get current bf
                BF& bf = vBF[curr_bf];
                
                //	set node id == scv this bf belongs to
                if (m_roid != ROID_PYRAMID || side != 0){
                    bf.nodeId = rRefElem.id(dim-1, side, 0, co);
                }
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
                traits::NormalOnBF(bf.Normal, bf.vGloPos, m_vvGloMid[0]);
                    
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
                    
            } // end loop of corners of side
                
        } // end loop sides of element
     
    }
    UG_CATCH_THROW("DimFV1FTGeometry: update_inner_boundary() failed.");
    
    
}

// compare implementation of 'DimFV1Geometry::update_boundary_faces()'
template <int TDim, int TWorldDim, class TInterfaceHandler>
void DimFV1FTGeometry<TDim, TWorldDim,  TInterfaceHandler>::
update_inner_boundary_faces_for2()
{
    for ( size_t i = 0; i < 4; ++i )
        //UG_LOG("m_spInterfaceHandler->quadriCorners()[" << i << "]: " << m_spInterfaceHandler->quadriCorners()[i] << "\n");
            
    /////////////////////////////////////////////////////////////////////////////
    // get general data
    /////////////////////////////////////////////////////////////////////////////
    try{
        const DimReferenceElement<dim>& rRefElem
        = ReferenceElementProvider::get<dim>(ROID_QUADRILATERAL);
                
        DimReferenceMapping<dim, worldDim>& rMapping = ReferenceMappingProvider::get<dim, worldDim>(ROID_QUADRILATERAL);
        rMapping.update(m_spInterfaceHandler->quadriCorners());
                
        const LocalShapeFunctionSet<dim>& TrialSpace =
        LocalFiniteElementProvider::get<dim>(ROID_QUADRILATERAL, LFEID(LFEID::LAGRANGE, dim, 1));
                
        MathVector<dim> buffer_vvLocMid[dim+1][maxMid];
        // 	set corners of element as local centers of nodes
        for(size_t i = 0; i < rRefElem.num(0); ++i)
        {
            buffer_vvLocMid[0][i] = rRefElem.corner(i);
            UG_LOG("**** buffer_vvLocMid[0][" << i << "]: " << buffer_vvLocMid[0][i] << "\n");
            
        }
        //	compute local midpoints
        ComputeMidPoints<dim, DimReferenceElement<dim>, maxMid>
        (rRefElem, buffer_vvLocMid[0], buffer_vvLocMid);
        
        MathVector<worldDim> buffer_vvGloMid[dim+1][maxMid];
        // 	remember global position of nodes
        for(size_t i = 0; i < rRefElem.num(0); ++i)
        {
            buffer_vvGloMid[0][i] = m_spInterfaceHandler->quadriCorners()[i];
            UG_LOG("**** buffer_vvGloMid[0][" << i << "]: " << buffer_vvGloMid[0][i] << "\n");
        }
        //	compute local midpoints
        ComputeMidPoints<worldDim, DimReferenceElement<dim>, maxMid>(rRefElem, buffer_vvGloMid[0], buffer_vvGloMid);
                
                
        /////////////////////////////////////////////////////////////////////////////
        // collect boudary faces
        /////////////////////////////////////////////////////////////////////////////
                
        // get number of sides of element
        size_t numSides = 0;
        numSides = rRefElem.num(dim-1);
        
        //	current number of bf
        size_t curr_bf = 0;
        
        std::vector<BF>& vBF = m_spInterfaceHandler->get_boundary_faces();
        
        vBF.clear();
        //	loop sides of element
        for(size_t side = 0; side < numSides; ++side)
        {
            // side is no boundary face => continue
            if ( side == 1 || side == 3 )
                continue;
                    
            //	number of corners of side (special case bottom side pyramid)
            const int coOfSide = (ROID_QUADRILATERAL != ROID_PYRAMID || side != 0)
            ? rRefElem.num(dim-1, side, 0) : rRefElem.num(dim-1, side, 0) + 2;
                    
            //	resize vector
            vBF.resize(curr_bf + coOfSide);
                    
            //	loop corners
            for(int co = 0; co < coOfSide; ++co)
            {
                //	get current bf
                BF& bf = vBF[curr_bf];
                        
                // set nodeID to ID of corner within the Triangle used for assembling the fluid node equation:
                bf.nodeId = curr_bf; //m_spInterfaceHandler->m_vQuadriOrigID[curr_bf]; //rRefElem.id(dim-1, side, 0, co);
                        
                //	Compute MidID for BF
                        ComputeBFMidID(rRefElem, side, bf.vMidID, co);
            
                // 	copy corners of bf
                CopyCornerByMidID<dim, maxMid>(bf.vLocPos, bf.vMidID, buffer_vvLocMid, BF::numCo);
                CopyCornerByMidID<worldDim, maxMid>(bf.vGloPos, bf.vMidID, buffer_vvGloMid, BF::numCo);
                
                // 	integration point
                AveragePositions(bf.localIP, bf.vLocPos, BF::numCo);
                AveragePositions(bf.globalIP, bf.vGloPos, BF::numCo);
                        
                // 	normal on scvf
                traits::NormalOnBF(bf.Normal, bf.vGloPos, buffer_vvGloMid[0]);
                        
                UG_LOG("BF::numCo " << BF::numCo << "\n");
                UG_LOG("co = " << co << ": bf.vGloPos[0]: " << bf.vGloPos[0] << "\n");
                UG_LOG("co = " << co << ": bf.vGloPos[1]: " << bf.vGloPos[1] << "\n");
                UG_LOG("bf.Normal:"  << bf.Normal <<  "\n");
                UG_LOG("bf.nodeID:"  << bf.nodeId <<  "\n");
                        
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
                
            } // end loop of corners of side
            
        } // end loop sides of element
        
    }
    UG_CATCH_THROW("DimFV1FTGeometry: update_inner_boundary() failed.");
        
}

    
template <int TDim, int TWorldDim, class TInterfaceHandler>
void DimFV1FTGeometry<TDim, TWorldDim,  TInterfaceHandler>::
update_boundary_faces(GridObject* pElem, const MathVector<worldDim>* vCornerCoords, const ISubsetHandler* ish)
{
//	get grid
	Grid& grid = *(ish->grid());

//	vector of subset indices of side
	std::vector<int> vSubsetIndex;

//	get subset indices for sides (i.e. edge in 2d, faces in 3d)
	if(dim == 1) {
		std::vector<Vertex*> vVertex;
		CollectVertices(vVertex, grid, pElem);
		vSubsetIndex.resize(vVertex.size());
		for(size_t i = 0; i < vVertex.size(); ++i)
			vSubsetIndex[i] = ish->get_subset_index(vVertex[i]);
	}
	if(dim == 2) {
		std::vector<Edge*> vEdges;
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
		LocalFiniteElementProvider::get<dim>(m_roid, LFEID(LFEID::LAGRANGE, dim, 1));

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

		//	number of corners of side
			const int coOfSide = rRefElem.num(dim-1, side, 0);

		//	resize vector
			vBF.resize(curr_bf + coOfSide);

		//	loop corners
			for(int co = 0; co < coOfSide; ++co)
			{
			//	get current bf
				BF& bf = vBF[curr_bf];

			//	set node id == scv this bf belongs to
				bf.nodeId = rRefElem.id(dim-1, side, 0, co);

			//	Compute MidID for BF
				ComputeBFMidID(rRefElem, side, bf.vMidID, co);

			// 	copy corners of bf
				CopyCornerByMidID<dim, maxMid>(bf.vLocPos, bf.vMidID, m_vvLocMid, BF::numCo);
				CopyCornerByMidID<worldDim, maxMid>(bf.vGloPos, bf.vMidID, m_vvGloMid, BF::numCo);

			// 	integration point
				AveragePositions(bf.localIP, bf.vLocPos, BF::numCo);
				AveragePositions(bf.globalIP, bf.vGloPos, BF::numCo);

			// 	normal on scvf
				traits::NormalOnBF(bf.Normal, bf.vGloPos, m_vvGloMid[0]);

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
	UG_CATCH_THROW("DimFV1FTGeometry: update failed.");
}

//////////////////////
// FV1FTGeometry

template class FV1FTGeometry<RegularEdge, 1>;
template class FV1FTGeometry<RegularEdge, 2>;
template class FV1FTGeometry<RegularEdge, 3>;

template class FV1FTGeometry<Triangle, 2>;
template class FV1FTGeometry<Triangle, 3>;

template class FV1FTGeometry<Quadrilateral, 2>;
template class FV1FTGeometry<Quadrilateral, 3>;

template class FV1FTGeometry<Tetrahedron, 3>;
template class FV1FTGeometry<Prism, 3>;
template class FV1FTGeometry<Pyramid, 3>;
template class FV1FTGeometry<Hexahedron, 3>;
template class FV1FTGeometry<Octahedron,3>;

//////////////////////
// DimFV1FTGeometry
template class DimFV1FTGeometry<1, 1, InterfaceHandlerLocalParticle<1> >;
    
template class DimFV1FTGeometry<2, 2, InterfaceHandlerLocalParticle<2> >;
    
template class DimFV1FTGeometry<3, 3, InterfaceHandlerLocalParticle<3> >;
    
    
} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DISC_HELPER__FV1FT_GEOMETRY_IMPL_ */

/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Martin Scherer, Rebecca Wittum, Sebastian Reiter
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

#include <cmath>
#include "tkd_info.h"
#include "common/math/misc/math_constants.h"

namespace ug {

/** Each line corresponds to one element. First entry of each line corresponds
 * to the number of vertices, the other entries to the vertex index according
 * to TKDInfo.*/
static const int INNER_TKD_ELEMENT_INDICES[] = {
					6, 24, 1, 0, 25, 19, 18,
					6, 24, 2, 1, 25, 20, 19,
					6, 24, 3, 2, 25, 21, 20,
					6, 24, 4, 3, 25, 22, 21,
					6, 24, 5, 4, 25, 23, 22,
					6, 24, 0, 5, 25, 18, 23,

					6, 2, 14, 20, 3, 15, 21,
					6, 3, 9, 21, 4, 10, 22,
					6, 4, 16, 22, 5, 17, 23,
					6, 5, 11, 23, 0, 6, 18,
					6, 0, 12, 18, 1, 13, 19,
					6, 1, 7, 19, 2, 8, 20,

					4, 3, 21, 9, 15,
					4, 4, 22, 10, 16,
					4, 5, 23, 11, 17,
					4, 0, 18, 6, 12,
					4, 1, 19, 7, 13,
					4, 2, 20, 8, 14};


/** Each line corresponds to one element. First entry of each line corresponds
 * to the number of vertices, the other entries to the vertex index according
 * to TKDInfo.*/
static const int OUTER_TKD_ELEMENT_INDICES[] = {
					//top
					6, 24, 1, 0, 26 + 24, 26 + 1, 26 + 0,
					6, 24, 2, 1, 26 + 24, 26 + 2, 26 + 1,
					6, 24, 3, 2, 26 + 24, 26 + 3, 26 + 2,
					6, 24, 4, 3, 26 + 24, 26 + 4, 26 + 3,
					6, 24, 5, 4, 26 + 24, 26 + 5, 26 + 4,
					6, 24, 0, 5, 26 + 24, 26 + 0, 26 + 5,

					//bottom
					6, 25, 18, 19, 26 + 25, 26 + 18, 26 + 19,
					6, 25, 19, 20, 26 + 25, 26 + 19, 26 + 20,
					6, 25, 20, 21, 26 + 25, 26 + 20, 26 + 21,
					6, 25, 21, 22, 26 + 25, 26 + 21, 26 + 22,
					6, 25, 22, 23, 26 + 25, 26 + 22, 26 + 23,
					6, 25, 23, 18, 26 + 25, 26 + 23, 26 + 18,

					//outer edges
					6, 0, 12, 6, 26 + 0, 26 + 12, 26 + 6,
					6, 1, 7, 13, 26 + 1, 26 + 7, 26 + 13,
					6, 2, 14, 8, 26 + 2, 26 + 14, 26 + 8,
					6, 3, 9, 15, 26 + 3, 26 + 9, 26 + 15,
					6, 4, 16, 10, 26 + 4, 26 + 16, 26 + 10,
					6, 5, 11, 17, 26 + 5, 26 + 11, 26 + 17,

					6, 18, 6, 12, 26 + 18, 26 + 6, 26 + 12,
					6, 19, 13, 7, 26 + 19, 26 + 13, 26 + 7,
					6, 20, 8, 14, 26 + 20, 26 + 8, 26 + 14,
					6, 21, 15, 9, 26 + 21, 26 + 15, 26 + 9,
					6, 22, 10, 16, 26 + 22, 26 + 10, 26 + 16,
					6, 23, 17, 11, 26 + 23, 26 + 17, 26 + 11,

					// hexahedra for quads
					8, 1, 2, 8, 7, 26 + 1, 26 + 2, 26 + 8, 26 + 7,
					8, 3, 4, 10, 9, 26 + 3, 26 + 4, 26 + 10, 26 + 9,
					8, 5, 0, 6, 11, 26 + 5, 26 + 0, 26 + 6, 26 + 11,

					8, 0, 1, 13, 12, 26 + 0, 26 + 1, 26 + 13, 26 + 12,
					8, 2, 3, 15, 14, 26 + 2, 26 + 3, 26 + 15, 26 + 14,
					8, 4, 5, 17, 16, 26 + 4, 26 + 5, 26 + 17, 26 + 16,

					8, 7, 8, 20, 19, 26 + 7, 26 + 8, 26 + 20, 26 + 19,
					8, 9, 10, 22, 21, 26 + 9, 26 + 10, 26 + 22, 26 + 21,
					8, 11, 6, 18, 23, 26 + 11, 26 + 6, 26 + 18, 26 + 23,

					8, 12, 13, 19, 18, 26 + 12, 26 + 13, 26 + 19, 26 + 18,
					8, 14, 15, 21, 20, 26 + 14, 26 + 15, 26 + 21, 26 + 20,
					8, 16, 17, 23, 22, 26 + 16, 26 + 17, 26 + 23, 26 + 22};



TKDInfo::
TKDInfo (number a, number w, number h, number d)
{
	number beta = 1.0 / 2.0 * PI + acos(h / sqrt(h*h+3.0*(w-2.0*a)*(w-2.0*a)));
	number gamma = 1.0 / 2.0 * PI + acos(2.0*h / sqrt(4.0*h*h+3.0*(2.0*w-a)*(2.0*w-a)));
	number m1 = d / (2.0 * tan(beta / 2.0));
	number m2 = d / (2.0 * tan(gamma / 2.0));

	number a2 = sqrt(3) / 3.0 * (sqrt(3) * a + m1 + m2);
	number h2 = h+d;
	number w2 = h2 * (w-2*a)/h + 2.0*a2;

	m_coords.resize(NUM_COORDS);
	init_coords (&m_coords[0], a, w, h);
	init_coords (&m_coords[NUM_INNER_COORDS], a2, w2, h2);
}


void TKDInfo::
init_coords (	vector3* coordsOut,
				number a,
				number w,
				number h)
{
	//Coordinates of TOP points
	coordsOut[0] = vector3(-a/2, sqrt(3)*a/2, h/2);
	coordsOut[1] = vector3(a/2, sqrt(3)*a/2, h/2);
	coordsOut[2] = vector3(a, 0, h/2);
	coordsOut[3] = vector3(a/2, -sqrt(3)*a/2, h/2);
	coordsOut[4] = vector3(-a/2, -sqrt(3)*a/2, h/2);
	coordsOut[5] = vector3(-a, 0, h/2);

	//Coordinates of point on 2/3 of TKD h
	coordsOut[6] = vector3(-(w-a)/2, sqrt(3)/6*(w+a), h/6);
	coordsOut[7] = vector3((w-a)/2, sqrt(3)/6*(w+a), h/6);
	coordsOut[8] = vector3(w/2, sqrt(3)/6*(w-2*a), h/6);
	coordsOut[9] = vector3(a/2, -sqrt(3)/3*(w-a/2), h/6);
	coordsOut[10] = vector3(-a/2, -sqrt(3)/3*(w-a/2), h/6);
	coordsOut[11] = vector3(-w/2, sqrt(3)/6*(w-2*a), h/6);

	//Coordinates of point on 1/3 of TKD h
	coordsOut[12] = vector3(-a/2, +sqrt(3)/3*(w-a/2), -h/6);
	coordsOut[13] = vector3(+a/2, +sqrt(3)/3*(w-a/2), -h/6);
	coordsOut[14] = vector3(w/2, -sqrt(3)/6*(w-2*a), -h/6);
	coordsOut[15] = vector3((w-a)/2, -sqrt(3)/6*(w+a), -h/6);
	coordsOut[16] = vector3(-(w-a)/2, -sqrt(3)/6*(w+a), -h/6);
	coordsOut[17] = vector3(-w/2, -sqrt(3)/6*(w-2*a), -h/6);

	//Coordinates of BOT points
	coordsOut[18] = vector3(-a/2, sqrt(3)*a/2, -h/2);
	coordsOut[19] = vector3(a/2, sqrt(3)*a/2, -h/2);
	coordsOut[20] = vector3(a, 0, -h/2);
	coordsOut[21] = vector3(a/2, -sqrt(3)*a/2, -h/2);
	coordsOut[22] = vector3(-a/2, -sqrt(3)*a/2, -h/2);
	coordsOut[23] = vector3(-a, 0, -h/2);

	//Coordinates of centre points of TOP and BOT
	coordsOut[24] = vector3(0, 0, h/2);
	coordsOut[25] = vector3(0, 0, -h/2);
}

const int* TKDInfo::
inner_element_indices () const
{
	return INNER_TKD_ELEMENT_INDICES;
}


const int* TKDInfo::
outer_element_indices () const
{
	return OUTER_TKD_ELEMENT_INDICES;
}


}//	end of namespace

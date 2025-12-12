/*
 * Copyright (c) 2017:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#include "hexahedron_util.h"

#include "../grid_objects/prism_rules.h"
#include "geom_obj_util/face_util.h"

namespace ug {

Hexahedron* CreateHexahedronFromPrisms(Grid& grid, Prism* p0, Prism* p1)
{
//	find the quad of p0 which is contained in p1
	const int* quadInds = prism_rules::QUADS;

	int iq[2] = {-1, -1};
	for(int iquad = 0; iquad < 3; ++iquad){
		Face* q0 = grid.get_face(p0, quadInds[iquad]);
		iq[1] = GetFaceIndex(p1, q0);
		if(iq[1] != -1){
			iq[0] = quadInds[iquad];
			break;
		}
	}

	if(iq[0] == -1)
		return nullptr;

//	we mark for both prisms which vertices are shared
	bool isShared[2][6]	= {	{false, false, false, false, false, false},
							{false, false, false, false, false, false}};

	for(int iprism = 0; iprism < 2; ++iprism){
		const int* vrts = prism_rules::FACE_VRT_INDS[iq[iprism]];
		for(int ivrt = 0; ivrt < 4; ++ivrt)
			isShared[iprism][vrts[ivrt]] = true;
	}

//	access the vertices of both prisms
	Prism::ConstVertexArray vrts0 = p0->vertices();
	Prism::ConstVertexArray vrts1 = p1->vertices();

//	mark the two shared vertices in the bottom face of p0
	grid.begin_marking();
	int ubv[2] = {-1, -1};	//	unmarked bottom vertex index
	for(int ivrt = 0; ivrt < 3; ++ivrt){
		if(isShared[0][ivrt])
			grid.mark(vrts0[ivrt]);
		else
			ubv[0] = ivrt;
	}

//	find the matching bottom face in p1
	for(int itri = 0; itri < 2; ++itri){
		bool gotMarked = false;
		int nonMarked = -1;
		for(int ivrt = 0; ivrt < 3; ++ivrt){
			int index = prism_rules::FACE_VRT_INDS[prism_rules::TRIS[itri]][ivrt];
			if(grid.is_marked(vrts1[index]))
				gotMarked = true;
			else
				nonMarked = index;
		}

		if(gotMarked){
			ubv[1] = nonMarked;
			break;
		}
	}

	grid.end_marking();

	if(ubv[0] == -1 || ubv[1] == -1)
		return nullptr;

	const int base[4] = {ubv[0], (ubv[0] + 1) % 3, ubv[1], (ubv[0] + 2) % 3};

	return *grid.create<Hexahedron>(
			HexahedronDescriptor(
				vrts0[base[0]], vrts0[base[1]], vrts1[base[2]], vrts0[base[3]],
				vrts0[base[0] + 3], vrts0[base[1] + 3], vrts1[(base[2] + 3) % 6], vrts0[base[3] + 3]),
			p0);
}

}//	end of namespace

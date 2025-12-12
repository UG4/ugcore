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

#include "tkd_util.h"

#include "lib_grid/grid_objects/grid_objects.h"

using namespace std;

namespace ug {

void CreateTKDVertices (const TKDInfo& tkdInfo,
                        Grid& g,
                        APosition3& aPos,
                        Vertex** vrtsOut)
{
	Grid::VertexAttachmentAccessor aaPos(g, aPos);

	const vector3* coords = tkdInfo.coords();
	for(size_t i = 0; i < tkdInfo.num_coords(); ++i){
		Vertex* vrt = *g.create<RegularVertex>();
		aaPos[vrt] = coords[i];
		if(vrtsOut)
			vrtsOut[i] = vrt;
	}
}


/** \param volsOut	(optional) Array of size numElements.*/
static
void CreateVolumesFromElementIndexList (
            Grid& g,
			const int* elemIndexList,
			int numElements,
			Vertex* const* vrts,
			Volume** volsOut = nullptr)
{
	VolumeDescriptor vd;

	const int* inds = elemIndexList;

	for(int i = 0; i < numElements; ++i){
		const int elemId = *inds++;
		const size_t num = (size_t)elemId;
		vd.set_num_vertices(num);
		for(size_t j = 0; j < num; ++j){
			vd.set_vertex(j, vrts[*inds++]);
		}

		Volume* vol = nullptr;
		switch(elemId){
			case 4:	vol = *g.create<Tetrahedron>(TetrahedronDescriptor(vd)); break;
			case 5:	vol = *g.create<Pyramid>(PyramidDescriptor(vd)); break;
			case 6:	vol = *g.create<Prism>(PrismDescriptor(vd)); break;
			case 8:	vol = *g.create<Hexahedron>(HexahedronDescriptor(vd));	break;
			default: UG_THROW("Unsupported element index: " << elemId);
		}

		if(volsOut)
			volsOut[i] = vol;
	}
}


void CreateTKD (const TKDInfo& tkdInfo,
                Grid& g,
				APosition3& aPos,
				Volume** volsOut)
{
	Vertex* vrts[TKDInfo::NUM_COORDS];
	CreateTKDVertices (tkdInfo, g, aPos, vrts);
	CreateVolumesFromElementIndexList (g, tkdInfo.inner_element_indices(),
	                                   TKDInfo::NUM_INNER_ELEMENTS,
	                                   vrts,
	                                   volsOut);
}


void CreateTKDWithOuterLayer (const TKDInfo& tkdInfo,
	                          Grid& g,
	                          APosition3& aPos,
							  Volume** volsOut)
{
	Vertex* vrts[TKDInfo::NUM_COORDS];
	CreateTKDVertices (tkdInfo, g, aPos, vrts);
	CreateVolumesFromElementIndexList (g, tkdInfo.inner_element_indices(),
	                                   TKDInfo::NUM_INNER_ELEMENTS,
	                                   vrts,
	                                   volsOut);
	CreateVolumesFromElementIndexList (g, tkdInfo.outer_element_indices(),
	                                   TKDInfo::NUM_OUTER_ELEMENTS,
	                                   vrts,
	                                   volsOut ? volsOut + TKDInfo::NUM_INNER_ELEMENTS : 0);
}

}//	end of namespace

/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Authors: Sebastian Reiter, Nicolas Tessore
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

#include "file_io_ng.h"

#include <vector>

#include "lib_grid/lg_base.h"

extern "C" {
#include <lib_grid/file_io/externals/include/ng/ng.h>
}

using namespace std;

enum ng_face_type
{
	ng_triangle = 3,
	ng_quadrilateral = 4,
};

enum ng_volume_type
{
	ng_tetrahedra = 4,
	ng_pyramid = 5,
	ng_prism = 6,
	ng_hexahedra = 8
};

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	ImportGridFromNG
bool ImportGridFromNG(Grid& grid,
                      const char* filename,
                      AVector3& aPos,
                      ISubsetHandler* pSubdomainHandler)
{
	// create ng object
	ng* n = ng_new();

	// read ng file
	ng_info* ninfo = ng_info_new();
	if(ng_read(filename, n, ninfo))
	{
		ng_info_delete(ninfo);
		ng_delete(n);
		n = ng_new();
		n->dim = 2;
		ninfo = ng_info_new();
		
	//TODO: 2d parsing fails, since it is not fully supported.
	//		There are problems in the element description (awaits F).
		if(ng_read(filename, n, ninfo)){
			LOG("WARNING in ImportGridFromNG: " << ninfo->err_msg << endl);
			ng_info_delete(ninfo);
			ng_delete(n);
			return false;
		}
	}
	ng_info_delete(ninfo);

	//	set up vertex attachment
	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);

	//	set up vertex attachment accessor
	Grid::VertexAttachmentAccessor aaPosition(grid, aPos);

	//	read nodes and store them in an array for index access
	vector<RegularVertex*>	vVertices;
	vVertices.reserve(n->num_bnodes + n->num_inodes);

	// read boundary nodes
	for(int i = 0; i < n->num_bnodes; ++i)
	{
		// get boundary node
		ng_bnode* node = &n->bnodes[i];

		// create and store vertex
		RegularVertex* vert = *grid.create<RegularVertex>();
		vVertices.push_back(vert);

		// set vertex coordinates
		aaPosition[vert] = vector3(
			node->coords[0],
			node->coords[1],
			node->coords[2]
		);
	}

	// read interior nodes
	for(int i = 0; i < n->num_inodes; ++i)
	{
		// get interior node
		ng_inode* node = &n->inodes[i];

		// create and store vertex
		RegularVertex* vert = *grid.create<RegularVertex>();
		vVertices.push_back(vert);

		// set vertex coordinates
		aaPosition[vert] = vector3(
			node->coords[0],
			node->coords[1],
			node->coords[2]
		);
	}

//	if we're in 2d, set all z-coords to 0
	if(n->dim == 2){
		for(size_t i = 0; i < vVertices.size(); ++i)
			aaPosition[vVertices[i]].z() = 0;
	}

//	create the elements
	if(n->dim == 2){
		//	read faces
		for(int i = 0; i < n->num_elements; ++i)
		{
			ng_element* elem = &n->elements[i];

			Face* face = nullptr;

			// create face
			switch(elem->num_nodes)
			{
				case ng_face_type::ng_triangle:
					face = *grid.create<Triangle>(TriangleDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]]));
					break;

				case ng_face_type::ng_quadrilateral:
					face = *grid.create<Quadrilateral>(QuadrilateralDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]],
							vVertices[elem->nodes[3]]));
					break;
				default:
					LOG("WARNING in ImportGridFromNG: Face type not implemented!" << endl);
					break;
			}

			// add face to subset
			if(face != nullptr && pSubdomainHandler != nullptr)
				pSubdomainHandler->assign_subset(face, elem->subdomain - 1);
		}
	}
	else{
		//	read volumes
		for(int i = 0; i < n->num_elements; ++i)
		{
			ng_element* elem = &n->elements[i];

			Volume* vol = nullptr;

			// create volume
			switch(elem->num_nodes)
			{
				case ng_volume_type::ng_tetrahedra:
					vol = *grid.create<Tetrahedron>(TetrahedronDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]],
							vVertices[elem->nodes[3]]
					));
					break;

				case ng_volume_type::ng_pyramid:
					vol = *grid.create<Pyramid>(PyramidDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]],
							vVertices[elem->nodes[3]],
							vVertices[elem->nodes[4]]));
					break;
				case ng_volume_type::ng_prism:
					vol = *grid.create<Prism>(PrismDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]],
							vVertices[elem->nodes[3]],
							vVertices[elem->nodes[4]],
							vVertices[elem->nodes[5]]));
					break;
				case ng_volume_type::ng_hexahedra:
					vol = *grid.create<Hexahedron>(HexahedronDescriptor(
							vVertices[elem->nodes[0]],
							vVertices[elem->nodes[1]],
							vVertices[elem->nodes[2]],
							vVertices[elem->nodes[3]],
							vVertices[elem->nodes[4]],
							vVertices[elem->nodes[5]],
							vVertices[elem->nodes[6]],
							vVertices[elem->nodes[7]]));
					break;
				default:
					LOG("WARNING in ImportGridFromNG: Volume type not implemented!" << endl);
					break;
			}

			// add volume to subset
			if(vol != nullptr && pSubdomainHandler != nullptr)
				pSubdomainHandler->assign_subset(vol, elem->subdomain - 1);
		}
	}
	// done importing!

	// delete ng object
	ng_delete(n);

	return true;
}

}//	end of namespace

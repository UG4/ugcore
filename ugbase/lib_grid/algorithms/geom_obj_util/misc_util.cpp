/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
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

//ø #include "common/util/vec_for_each.h"
#include "misc_util.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	EraseConnectingElements
///	erases all elements that connect v1 and v2
void EraseConnectingElements(Grid& grid, Vertex* v1, Vertex* v2)
{
//	check volumes for direct connection
	if(grid.num_volumes() > 0)
	{
	//	iterate through associated volumes
		std::vector<Volume*> vols;
		CollectVolumes(vols, grid, v1);
		for(size_t _vfeI = 0; _vfeI < vols.size(); ++_vfeI){
			Volume* v = vols[_vfeI];
			uint numVrts = v->num_vertices();
			bool gotOne = false;
			for(uint i = 0; i < numVrts; ++i)
			{
				if(v->vertex(i) == v2)
				{
					gotOne = true;
					break;
				}
			}

			if(gotOne == true)
				grid.erase(v);
		}
	}
	
//	check faces for direct connection.
	if(grid.num_faces() > 0)
	{
	//	iterate through associated faces
		std::vector<Face*> faces;
		CollectFaces(faces, grid, v1);
		for(size_t _vfeI = 0; _vfeI < faces.size(); ++_vfeI){
			Face* f = faces[_vfeI];
			uint numVrts = f->num_vertices();
			bool gotOne = false;
			for(uint i = 0; i < numVrts; ++i)
			{
				if(f->vertex(i) == v2)
				{
					gotOne = true;
					break;
				}
			}

			if(gotOne == true)
				grid.erase(f);
		}
	}

//	check edges
	if(grid.num_edges() > 0)
	{
	//	iterate through associated edges
		std::vector<Edge*> edges;
		CollectEdges(edges, grid, v1);
		for(size_t _vfeI = 0; _vfeI < edges.size(); ++_vfeI){
			Edge* e = edges[_vfeI];
		//	if e contains v2 we have to remove it.
			if((e->vertex(0) == v2) || (e->vertex(1) == v2))
			{
				grid.erase(e);
			}
		}
	}
}

}//	end of namespace

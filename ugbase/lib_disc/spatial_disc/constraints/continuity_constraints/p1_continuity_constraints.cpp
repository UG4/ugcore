/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "p1_continuity_constraints.h"

//#include "lib_grid/algorithms/geom_obj_util/edge_util.h"
#include "lib_grid/algorithms/debug_util.h"

namespace ug {


/// returns the vertices of the object constraining a hanging vertex
void CollectConstraining(std::vector<Vertex*>& vConstrainingVrt,
						 const Grid& grid,
                         ConstrainedVertex* hgVrt,
                         bool bClearContainer)
{
//	clear container
	if(bClearContainer) vConstrainingVrt.clear();

//	switch constraining parent
	switch(hgVrt->get_parent_base_object_id())
	{
	case EDGE:
	{
	//	in a parallel environment, the parent may be missing...
		GridObject* constrainingObject = hgVrt->get_constraining_object();
		if(constrainingObject){
		//	cast to constraining edge
			auto constrainingEdge = dynamic_cast<ConstrainingEdge*>(constrainingObject);

		//	check that edge is correct
			if(constrainingEdge == nullptr)
				UG_THROW("Parent element should be "
							"constraining edge, but is not.");

		//	get constraining vertices
			for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
			{
			//	get constrained edge
				auto constrainedEdge = dynamic_cast<ConstrainedEdge*>( constrainingEdge->constrained_edge(i_cde));

			//	check
				if(constrainedEdge == nullptr)
					UG_THROW("Child element should be "
								"constrained edge, but is not.");

			//	get non-hanging vertex
				Vertex* vrt = GetConnectedVertex(constrainedEdge, hgVrt);

			//	push back in list of interpolation vertices
				vConstrainingVrt.push_back(vrt);
			}
		}
		else
		{
		//	we have to find the constraining vertices of hgVrt without having
		//	access to the parent element.
			Grid::edge_traits::secure_container edges;
		//todo: associated elements should support const!
			const_cast<Grid&>(grid).associated_elements(edges, hgVrt);
			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				Edge* e = edges[i_edge];
				if(e->is_constrained()){
					Vertex* conVrt = GetConnectedVertex(e, hgVrt);
					if(!conVrt->is_constrained()){
						vConstrainingVrt.push_back(conVrt);
					}
				}
			}
		}
	}
		break;
	case FACE:
	{
	//	in a parallel environment, the parent may be missing...
		GridObject* constrainingObject = hgVrt->get_constraining_object();
		if(constrainingObject){
		//	cast to constraining quadrilateral
			auto bigQuad = dynamic_cast<ConstrainingQuadrilateral*>(constrainingObject);

		//	check that quad is correct
			if(bigQuad == nullptr)
				UG_THROW("Parent element should be "
								"constraining quad, but is not.");

		//	get constraining vertices
		//	\todo: This is only valid for a surface grid!!!
			for(size_t i_cf=0; i_cf < bigQuad->num_constrained_faces(); ++i_cf)
			{
				Face* face = bigQuad->constrained_face(i_cf);

				Vertex* vrt = nullptr;
				size_t i_vrt = 0;
				for(i_vrt = 0; i_vrt < face->num_vertices(); ++i_vrt)
				{
					vrt = face->vertex(i_vrt);
					if(hgVrt != vrt && (!vrt->is_constrained()))
						break;
				}
				if(i_vrt == face->num_vertices())
					UG_THROW("ERROR: Vertex not detected.\n");

				vConstrainingVrt.push_back(vrt);
			}
		}
		else
		{
		//	we have to find the constraining vertices of hgVrt without having
		//	access to the parent element.
			Grid::face_traits::secure_container faces;
		//todo: associated elements should support const!
			const_cast<Grid&>(grid).associated_elements(faces, hgVrt);
			for(size_t i_face = 0; i_face < faces.size(); ++i_face){
				Face* f = faces[i_face];
				if(f->is_constrained()){
					Vertex* vrt = nullptr;
					size_t i_vrt = 0;
					for(i_vrt = 0; i_vrt < f->num_vertices(); ++i_vrt){
						vrt = f->vertex(i_vrt);
						if(hgVrt != vrt && (!vrt->is_constrained()))
							break;
					}

					if(i_vrt == f->num_vertices())
						UG_THROW("ERROR: Vertex not detected.\n");

					vConstrainingVrt.push_back(vrt);
				}
			}
		}
	}
		break;
	default:
		UG_THROW("Parent element of hanging vertex wrong: "
				<< hgVrt->get_parent_base_object_id()
				<< ". Element info: "
				<< ElementDebugInfo(grid, hgVrt));
		break;
	}
}

} // end namespace ug

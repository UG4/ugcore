/*
 * p1_continuity_constraints.cpp
 *
 *  Created on: 06.03.2012
 *      Author: andreasvogel
 */

#include "p1_continuity_constraints.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"

namespace ug{


/// returns the vertices of the object constraining a hanging vertex
void CollectConstraining(std::vector<VertexBase*>& vConstrainingVrt,
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
	//	in a parallel environment, the parent may be missing
		GeometricObject* constrainingObject = hgVrt->get_constraining_object();
		if(constrainingObject){
		//	cast to constraining edge
			ConstrainingEdge* constrainingEdge =
					dynamic_cast<ConstrainingEdge*>(constrainingObject);

		//	check that edge is correct
			if(constrainingEdge == NULL)
				UG_THROW("Parent element should be "
							"constraining edge, but is not.");

		//	get constraining vertices
			for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
			{
			//	get constrained edge
				ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
													constrainingEdge->constrained_edge(i_cde));

			//	check
				if(constrainedEdge == NULL)
					UG_THROW("Child element should be "
								"constrained edge, but is not.");

			//	get non-hanging vertex
				VertexBase* vrt = GetConnectedVertex(constrainedEdge, hgVrt);

			//	push back in list of interpolation vertices
				vConstrainingVrt.push_back(vrt);
			}
		}
		else{
		//	we have to find the constraining vertices of hgVrt without having
		//	access to the parent element.
			Grid::edge_traits::secure_container edges;
		//todo: associated elements should support const!
			const_cast<Grid&>(grid).associated_elements(edges, hgVrt);
			for(size_t i_edge = 0; i_edge < edges.size(); ++i_edge){
				EdgeBase* e = edges[i_edge];
				if(e->is_constrained()){
					VertexBase* conVrt = GetConnectedVertex(e, hgVrt);
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
		GeometricObject* constrainingObject = hgVrt->get_constraining_object();
		if(constrainingObject){
		//	cast to constraining quadrilateral
			ConstrainingQuadrilateral* bigQuad =
					dynamic_cast<ConstrainingQuadrilateral*>(constrainingObject);

		//	check that quad is correct
			if(bigQuad == NULL)
				UG_THROW("Parent element should be "
								"constraining quad, but is not.");

		//	get constraining vertices
		//	\todo: This is only valid for a surface grid!!!
			for(size_t i_cf=0; i_cf < bigQuad->num_constrained_faces(); ++i_cf)
			{
				Face* face = bigQuad->constrained_face(i_cf);

				VertexBase* vrt = NULL;
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
		else{
		//	we have to find the constraining vertices of hgVrt without having
		//	access to the parent element.
			Grid::face_traits::secure_container faces;
		//todo: associated elements should support const!
			const_cast<Grid&>(grid).associated_elements(faces, hgVrt);
			for(size_t i_face = 0; i_face < faces.size(); ++i_face){
				Face* f = faces[i_face];
				if(f->is_constrained()){
					VertexBase* vrt = NULL;
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
	default: UG_THROW("Parent element of hanging vertex wrong."); break;
	}
}

} // end namespace ug

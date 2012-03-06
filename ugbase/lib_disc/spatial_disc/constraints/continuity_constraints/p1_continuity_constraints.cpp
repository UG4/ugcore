/*
 * p1_continuity_constraints.cpp
 *
 *  Created on: 06.03.2012
 *      Author: andreasvogel
 */

#include "p1_continuity_constraints.h"

namespace ug{


/// returns the vertices of the object constraining a hanging vertex
void CollectConstraining(std::vector<VertexBase*>& vConstrainingVrt,
                         HangingVertex* hgVrt,
                         bool bClearContainer)
{
//	clear container
	if(bClearContainer) vConstrainingVrt.clear();

//	switch constraining parent
	switch(hgVrt->get_parent_base_object_type_id())
	{
	case EDGE:
	{
	//	cast to constraining edge
		ConstrainingEdge* constrainingEdge =
				dynamic_cast<ConstrainingEdge*>(hgVrt->get_parent());

	//	check that edge is correct
		if(constrainingEdge == NULL)
			UG_THROW_FATAL("Parent element should be "
						"constraining edge, but is not.");

	//	get constraining vertices
		for(size_t i_cde = 0; i_cde < constrainingEdge->num_constrained_edges(); ++i_cde)
		{
		//	get constrained edge
			ConstrainedEdge* constrainedEdge = dynamic_cast<ConstrainedEdge*>(
												constrainingEdge->constrained_edge(i_cde));

		//	check
			if(constrainedEdge == NULL)
				UG_THROW_FATAL("Child element should be "
							"constrained edge, but is not.");

		//	get non-hanging vertex
			VertexBase* vrt = GetConnectedVertex(constrainedEdge, hgVrt);

		//	push back in list of interpolation vertices
			vConstrainingVrt.push_back(vrt);
		}
	}
		break;
	case FACE:
	{
	//	cast to constraining quadrilateral
		ConstrainingQuadrilateral* bigQuad =
				dynamic_cast<ConstrainingQuadrilateral*>(hgVrt->get_parent());

	//	check that quad is correct
		if(bigQuad == NULL)
			UG_THROW_FATAL("Parent element should be "
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
				if(hgVrt != vrt && dynamic_cast<HangingVertex*>(vrt) == NULL)
					break;
			}
			if(i_vrt == face->num_vertices())
				UG_THROW_FATAL("ERROR: Vertex not detected.\n");

			vConstrainingVrt.push_back(vrt);
		}
	}
		break;
	default: UG_THROW_FATAL("Parent element of hang. vertex wrong.");
	}
}

} // end namespace ug

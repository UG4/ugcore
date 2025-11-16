/*
 * Copyright (c) 2016:  G-CSC, Goethe University Frankfurt
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

#include "quadrilateral_util.h"
#include "lib_grid/algorithms/geom_obj_util/edge_util.h"

namespace ug{

Quadrilateral* CreateQuadrilateral_NoRegistration(Grid& g, Face* tri1, Face* tri2)
{
	UG_COND_THROW(tri1->num_vertices() != 3 || tri2->num_vertices() != 3,
				  "Only 2 triangles can be converted to quadrilaterals, "
				  "but non-triangle specifed.");

	Edge* e = GetConnectingEdge(g, tri1, tri2);
	UG_COND_THROW(e == nullptr, "Only triangles sharing one edge may be converted "
				  "to a quadrilateral. This is not the case here.");

	Vertex* v1 = GetConnectedVertex(e, tri1);
	Vertex* v3 = GetConnectedVertex(e, tri2);

	UG_COND_THROW(v1 == nullptr || v3 == nullptr || v1 == v3,
				  "Invalid triangles specified. Only two triangles sharing "
				  "exactly 2 vertices can be converted to a triangle.");

	int i1 = GetVertexIndex(tri1, v1);
	Vertex* v2 = tri1->vertex((i1+1)%3);
	Vertex* v4 = tri1->vertex((i1+2)%3);

	return new Quadrilateral(v1, v2, v3, v4);
}

Quadrilateral* CreateQuadrilateral(Grid& g, Face* tri1, Face* tri2)
{
	Quadrilateral* q = CreateQuadrilateral_NoRegistration(g, tri1, tri2);
	g.register_element(q);
	return q;
}

Quadrilateral* ReplaceByQuadrilateral(Grid& g, Face* tri1, Face* tri2)
{
	Quadrilateral* q = CreateQuadrilateral_NoRegistration(g, tri1, tri2);
	Edge* e = GetConnectingEdge(g, tri1, tri2);
	g.erase(tri2);
	g.register_element(q, tri1);
	g.erase(tri1);
	g.erase(e);
	return q;
}
	
}//	end of namespace

/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "icosahedron.h"
#include "lib_grid/refinement/regular_refinement_new.h"
#include "lib_grid/refinement/projectors/sphere_projector.h"

namespace ug{

///	Creates an Icosahedron
void GenerateIcosahedron(Grid& grid, const vector3& center,
						 number radius, AVector3& aPos)
{
	const number A = 0.85065080835204;
	const number B = 0.525731112119134;

//	create the vertices
	const number coords[12][3] = {	{-B, A, 0}, {0, B, A}, {B, A, 0}, {0, B, -A},
									{-A, 0, B}, {A, 0, B}, {A, 0, -B}, {-A, 0, -B},
									{-B, -A, 0}, {0, -B, A}, {B, -A, 0}, {0, -B, -A}};

	const int inds[20][3] = {	{0, 1, 2}, {0, 2, 3},
								{3, 7, 0}, {7, 4, 0}, {0, 4, 1},
								{1, 5, 2}, {5, 6, 2}, {2, 6, 3},
								{4, 9, 1}, {1, 9, 5}, {7, 3, 11}, {3, 6, 11},
								{11, 8, 7}, {7, 8, 4}, {8, 9, 4},
								{9, 10, 5}, {10, 6, 5}, {10, 11, 6},
								{8, 10, 9}, {8, 11, 10}};

	Vertex* vrts[12];

	if(!grid.has_vertex_attachment(aPos))
		grid.attach_to_vertices(aPos);
	Grid::AttachmentAccessor<Vertex, AVector3> aaPos(grid, aPos);

	for(size_t i = 0; i < 12; ++i){
		vrts[i] = *grid.create<RegularVertex>();
		aaPos[vrts[i]] = vector3(coords[i][0], coords[i][1], coords[i][2]);
		VecScaleAdd(aaPos[vrts[i]], 1.0, center, radius, aaPos[vrts[i]]);
	}

//	construct triangles
	for(size_t i = 0; i < 20; ++i){
		grid.create<Triangle>(TriangleDescriptor(vrts[inds[i][0]], vrts[inds[i][1]],
												 vrts[inds[i][2]]));
	}
}

void GenerateIcosphere(Grid& grid, const vector3& center, number radius,
					   int numRefinements, AVector3& aPos, Selector* psel)
{
	Selector defSel;
	if(!psel){
		defSel.assign_grid(grid);
		psel = &defSel;
	}

	Selector& sel = *psel;

//	enable autoselection
	bool autoselectionEnabled = sel.autoselection_enabled();
	sel.enable_autoselection(true);

//	clear the selector
	sel.clear();

//	create an icosahedron. All elements of the sphere will be selected, since we
//	enabled autoselection
	GenerateIcosahedron(grid, center, radius, aPos);

//	perform refinement steps
//	we need a temporary int attachment for the edges
	AInt aInt;
	grid.attach_to_edges(aInt);

//	use a refinement callback to project the new vertices to the sphere
	SphereProjectorNew sphereProjecton(MakeGeometry3d(grid, aPos), center);

	for(int i = 0; i < numRefinements; ++i)
		RefineNew(grid, sel, aInt, &sphereProjecton);

//	remove attachments
	grid.detach_from_edges(aInt);

//	restore autoselection
	sel.enable_autoselection(autoselectionEnabled);
}

}// end of namespace

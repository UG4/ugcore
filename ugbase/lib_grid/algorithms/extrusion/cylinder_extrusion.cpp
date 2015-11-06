/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#include <vector>
#include "extrude.h"
#include "cylinder_extrusion.h"
#include "lib_grid/algorithms/remeshing/grid_adaption.h"
#include "lib_grid/algorithms/selection_util.h"

using namespace std;

namespace ug
{

//	this method actually performs the extrusion. It is only callable from
//	inside this file.
static bool ExtrudeCylinder(Grid& grid, SubsetHandler* sh, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd, int cylSubInd,
					Selector* pSel)
{
//	we'll use this selector if no selector was specified
	Selector locSel;
	if(!pSel){
		locSel.assign_grid(grid);
		pSel = &locSel;
	}
	Selector& sel = *pSel;
	
	if(!AdaptSurfaceGridToCylinder(sel, grid, vrt, direction, radius, rimSnapThreshold)){
		LOG("  WARNING: AdaptSurfaceGridToCylinder failed during ExtrudeCylinder.\n");
		return false;
	}

//	select boundary edges
	sel.clear<Edge>();
	SelectAreaBoundary(sel, sel.begin<Face>(), sel.end<Face>());

//	gather faces and edges for extrusion
	vector<Edge*> vEdges;
	vector<Face*> vFaces;

	for(EdgeIterator iter = sel.begin<Edge>();
		iter != sel.end<Edge>(); ++iter)
		vEdges.push_back(*iter);

	for(FaceIterator iter = sel.begin<Face>();
		iter != sel.end<Face>(); ++iter)
		vFaces.push_back(*iter);

//	everything that is extruded shall go into a separate subset
	bool bSInh = false;
	int defSubInd = -1;
	if(sh){
	//	assign faces that will be extruded into a separate subset
		sh->assign_subset(vFaces.begin(), vFaces.end(), bottomSubInd);
	//	the rest gets into the next subset
		bSInh = sh->subset_inheritance_enabled();
		defSubInd = sh->get_default_subset_index();
		sh->enable_subset_inheritance(false);
		sh->set_default_subset_index(cylSubInd);		
	}

	vector3 scaledDir;
	VecScale(scaledDir, direction, height);
	
	Extrude(grid, NULL, &vEdges, &vFaces, scaledDir, EO_CREATE_FACES, aPosition);

	if(sh){
	//	restore subset handler
		sh->enable_subset_inheritance(bSInh);
		sh->set_default_subset_index(defSubInd);
	}

	return true;
}

bool ExtrudeCylinder(Grid& grid, SubsetHandler& sh, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					int bottomSubInd, int cylSubInd,
					Selector* pSel)
{
	return ExtrudeCylinder(grid, &sh, vrt, direction, height, radius,
						 rimSnapThreshold, aaPos, bottomSubInd, cylSubInd, pSel);

}

bool ExtrudeCylinder(Grid& grid, Vertex* vrt,
					const vector3& direction, number height, number radius,
					number rimSnapThreshold,
					Grid::VertexAttachmentAccessor<APosition>& aaPos,
					Selector* pSel)
{
	return ExtrudeCylinder(grid, NULL, vrt, direction, height, radius,
							rimSnapThreshold, aaPos, -1, -1, pSel);
}
					
}//	end of namespace

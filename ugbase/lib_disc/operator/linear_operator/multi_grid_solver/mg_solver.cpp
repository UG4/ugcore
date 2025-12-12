/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
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

#include "mg_solver.h"

namespace ug {


////////////////////////////////////////////////////////////////////////////////
// SelectNonShadowsAdjacentToShadowsOnLevel
////////////////////////////////////////////////////////////////////////////////

void SelectNonShadowsAdjacentToShadowsOnLevel(BoolMarker& sel,
                                              const SurfaceView& surfView,
                                              int level)
{
	PROFILE_FUNC_GROUP("gmg");
//	vectors for associated elements
	std::vector<Vertex*> vAssVertex;
	std::vector<Edge*> vAssEdge;
	std::vector<Face*> vAssFace;
	std::vector<Volume*> vAssVolume;

//	get grid
	Grid& grid = *sel.grid();

//	get multigrid
	auto mg = dynamic_cast<MultiGrid*>(&grid);

//	check multigrid
	if(mg == nullptr)
		UG_THROW("SelectNonShadowsAdjacentToShadowsOnLevel: No "
					"Multigrid given, selection on level not possible.");

//	check level
	if(level >= static_cast<int>(mg->num_levels()) || level < 0)
		UG_THROW("SelectNonShadowsAdjacentToShadowsOnLevel: Requested "
						"level "<<level<<" does not exist in Multigrid.");

	const GridLevel gl(GridLevel::TOP, GridLevel::ViewType::SURFACE);

//	iterator type
	geometry_traits<Vertex>::const_iterator iter, iterEnd;

	iterEnd = mg->end<Vertex>(level);

//	loop all base elems
	for(iter = mg->begin<Vertex>(level); iter != iterEnd; ++iter)
	{
	//	get element
		Vertex* shadow = *iter;

	//	check if element is a shadow
		if(!surfView.is_shadowed(shadow)) continue;

	//	get adjacent elements
		CollectAssociated(vAssEdge, grid, shadow);

		vAssVertex.clear();
		for(size_t i = 0; i < vAssEdge.size(); ++i)
			CollectAssociated(vAssVertex, grid, vAssEdge[i], false);

		vAssEdge.clear();
		vAssFace.clear();
		vAssVolume.clear();
		for(size_t i = 0; i < vAssVertex.size(); ++i)
		{
			CollectAssociated(vAssEdge, grid, vAssVertex[i], false);
			CollectAssociated(vAssFace, grid, vAssVertex[i], false);
			CollectAssociated(vAssVolume, grid, vAssVertex[i], false);
		}

	//	select associated elements
		for(size_t i = 0; i < vAssVertex.size(); ++i)
			if(surfView.is_contained(vAssVertex[i], gl, SurfaceView::SURFACE))
				sel.mark(vAssVertex[i]);
		for(size_t i = 0; i < vAssEdge.size(); ++i)
			if(surfView.is_contained(vAssEdge[i], gl, SurfaceView::SURFACE))
				sel.mark(vAssEdge[i]);
		for(size_t i = 0; i < vAssFace.size(); ++i)
			if(surfView.is_contained(vAssFace[i], gl, SurfaceView::SURFACE))
				sel.mark(vAssFace[i]);
		for(size_t i = 0; i < vAssVolume.size(); ++i)
			if(surfView.is_contained(vAssVolume[i], gl, SurfaceView::SURFACE))
				sel.mark(vAssVolume[i]);
	}
}

} // end namespace ug

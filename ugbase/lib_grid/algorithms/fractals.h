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

#ifndef __H__UG__LIB_GRID__FRACTALS__
#define __H__UG__LIB_GRID__FRACTALS__

#include "lib_grid/algorithms/refinement/hanging_node_refiner_grid.h"
#include "lib_grid/algorithms/refinement/refinement_projectors/fractal_projector.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
///	Repeatedly refines a grid and moves new vertices along the normal by a given factor.
template <class TAPosition>
bool CreateFractal_NormalScale(Grid& grid, HangingNodeRefiner_Grid& href,
							   number scaleFac, size_t numIterations,
							   TAPosition& aPosVRT)
{
	if(!grid.has_vertex_attachment(aPosVRT)){
		UG_LOG("WARNING in CreateFractal_NormalScale: given position attachment is ");
		UG_LOG("not attached to the grids vertices. Aborting...\n");
		return false;
	}
		
//	store the old refinement callback of href
	IRefinementCallback* oldCallback = href.get_refinement_callback();
	
//	create the new one.
	FractalProjector refCallback(grid, scaleFac);
	href.set_refinement_callback(&refCallback);
	
//	iterate for the specified number of times
	for(size_t i = 0; i < numIterations; ++i){

		if(grid.num_volumes() > 0){
		//	iterate over all faces and mark them for refinement, if they are boundary faces.
			for(FaceIterator iter = grid.faces_begin();
				iter != grid.faces_end(); ++iter)
			{
				if(IsVolumeBoundaryFace(grid, *iter))
					href.mark(*iter);
			}
		}
		else if(grid.num_faces() > 0){
		//	markall faces
			href.mark(grid.faces_begin(), grid.faces_end());
		}
		else{
		//	mark all edges
			href.mark(grid.edges_begin(), grid.edges_end());
		}

	//	refine them
		href.refine();

	//	change the scalefac
		refCallback.set_scale_fac(-0.5 * refCallback.get_scale_fac());
		//refCallback.set_scale_fac(-0.5 * refCallback.get_scale_fac());
		//refCallback.set_scale_fac(refCallback.get_scale_fac() * refCallback.get_scale_fac());

	}

//	done. restore href
	href.set_refinement_callback(oldCallback);
	return true;
}

inline bool CreateFractal_NormalScale(Grid& grid,
									   HangingNodeRefiner_Grid& href,
									   number scaleFac, size_t numIterations)
{
	return CreateFractal_NormalScale(grid, href, scaleFac,
									 numIterations, aPosition);
}
							   
}

#endif

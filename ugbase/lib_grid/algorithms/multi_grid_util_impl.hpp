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

#ifndef __H__UG__MULTI_GRID_UTIL_IMPL__
#define __H__UG__MULTI_GRID_UTIL_IMPL__

#include "multi_grid_util.h"
#include "subset_util.h"

namespace ug {

template <typename TElem>
void CollectSurfaceViewElements(ISubsetHandler& surfaceViewOut,
								MultiGrid& mg,
								MultiGridSubsetHandler& mgsh,
								bool clearContainer)
{
	using ElemIter = typename geometry_traits<TElem>::iterator;

//	clear the target surfaceView
	if(clearContainer)
		surfaceViewOut.clear();
		
//	iterate through all levels of the mgsh
	for(size_t level = 0; level < mgsh.num_levels(); ++level){
//	iterate through all subsets on that level
		for(int si = 0; si < mgsh.num_subsets(); ++si){
		//	iterate through all elements in the subset on that level
			for(ElemIter iter = mgsh.begin<TElem>(si, level);
				iter != mgsh.end<TElem>(si, level); ++iter)
			{
				TElem* elem = *iter;
			//	check whether the element has children
				if(!mg.has_children(elem)){
				//	the element is a part of the surface-view. add it to the handler
					surfaceViewOut.assign_subset(elem, si);
				}
			}
		}
	}
}

template <typename TSurfaceView>
void CreateSurfaceView(TSurfaceView& surfaceViewOut,
						MultiGrid& mg,
						MultiGridSubsetHandler& mgsh)
{
//	This method clears the surfaceViewOut and assigns all objects of
//	which lie on the surface of the mg to the surface view.
	CollectSurfaceViewElements<Vertex>(surfaceViewOut, mg, mgsh, true);
	CollectSurfaceViewElements<Edge>(surfaceViewOut, mg, mgsh, false);
	CollectSurfaceViewElements<Face>(surfaceViewOut, mg, mgsh, false);
	CollectSurfaceViewElements<Volume>(surfaceViewOut, mg, mgsh, false);
	
//	assign associated elements of lower dimension to the surface view
	bool assignSidesOnly = true;
	if(mgsh.num<Volume>() > 0 && !mg.option_is_enabled(VolumeOptions::VOLOPT_AUTOGENERATE_FACES))
		assignSidesOnly = false;
	else if(mgsh.num<Volume>() > 0 && !mg.option_is_enabled(VolumeOptions::VOLOPT_AUTOGENERATE_EDGES))
		assignSidesOnly = false;
	else if(mgsh.num<Face>() > 0 && !mg.option_is_enabled(FaceOptions::FACEOPT_AUTOGENERATE_EDGES))
		assignSidesOnly = false;

	if(assignSidesOnly){
		AssignAssociatedSidesToSubsets<Volume>(surfaceViewOut, mgsh);
		AssignAssociatedSidesToSubsets<Face>(surfaceViewOut, mgsh);
		AssignAssociatedSidesToSubsets<Edge>(surfaceViewOut, mgsh);
	}
	else
	{
		UG_LOG("INFO in CreateSurfaceView: Performing AssignAssociatedLowerDimElemsToSubsets ");
		UG_LOG("for all elements (Small Performance drawback).\n");

		AssignAssociatedLowerDimElemsToSubsets<Volume>(surfaceViewOut, mgsh);
		AssignAssociatedLowerDimElemsToSubsets<Face>(surfaceViewOut, mgsh);
		AssignAssociatedLowerDimElemsToSubsets<Edge>(surfaceViewOut, mgsh);
	}
}

template <typename TElem>
bool IsSubSurfaceElement(MultiGrid& mg, TElem* e, bool checkSides)
{
	using TBaseElem = typename TElem::grid_base_object;

	size_t numChildren = mg.num_children<TBaseElem>(e);
	for(size_t i = 0; i < numChildren; ++i){
		TBaseElem* child = mg.get_child<TBaseElem>(e, i);
		assert(child);
		if(mg.num_children<TBaseElem>(child) != 0)
			return false;
	}

//	ok all children of the same base type are surface elements.
	if(checkSides){
	//	Now we have to check whether all sides are surface elements, too.
	//	Since vertices do not have sides, this recursion will terminate.
		for(size_t i = 0; i < e->num_sides(); ++i){
			if(!IsSubSurfaceElement(mg, mg.get_side(e, i)))
				return false;
		}
	}

	return true;
}

}//	end of namespace

#endif

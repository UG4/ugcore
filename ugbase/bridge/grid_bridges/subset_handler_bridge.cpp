/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
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

#include "grid_bridges.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/subset_handler.h"
#include "lib_grid/tools/surface_view.h"
using namespace std;

#ifdef UG_USE_PYBIND11
#include "bindings/pybind/ug_pybind.h"
#endif

namespace ug{

template <class TElem>
static void AssignSubsetsByLevel(SubsetHandler& sh, MultiGrid& mg)
{
	for(size_t lvl = 0; lvl < mg.num_levels(); ++lvl){
		std::stringstream ss;
		ss << "Level " << lvl;
		sh.subset_info(lvl).name = ss.str().c_str();

		typedef typename Grid::traits<TElem>::iterator TIter;
		for(TIter iter = mg.begin<TElem>(lvl); iter != mg.end<TElem>(lvl); ++iter){
			sh.assign_subset(*iter, lvl);
		}
	}
}

void AssignSubsetsByLevel(SubsetHandler& sh, MultiGrid& mg)
{
	AssignSubsetsByLevel<Vertex>(sh, mg);
	AssignSubsetsByLevel<Edge>(sh, mg);
	AssignSubsetsByLevel<Face>(sh, mg);
	AssignSubsetsByLevel<Volume>(sh, mg);

	AssignSubsetColors(sh);
	EraseEmptySubsets(sh);
}

namespace bridge{

template<typename TRegistry = Registry>
void RegisterGridBridge_SubsetHandler_(TRegistry& reg, string parentGroup)
{
	string grp = parentGroup;
	
//  ISubsetHandler
	reg.template add_class_<ISubsetHandler>("ISubsetHandler", grp)
		.add_method("assign_subset", static_cast<void (ISubsetHandler::*)(Vertex*, int)>(&ISubsetHandler::assign_subset))
		.add_method("assign_subset", static_cast<void (ISubsetHandler::*)(Edge*, int)>(&ISubsetHandler::assign_subset))
		.add_method("assign_subset", static_cast<void (ISubsetHandler::*)(Face*, int)>(&ISubsetHandler::assign_subset))
		.add_method("assign_subset", static_cast<void (ISubsetHandler::*)(Volume*, int)>(&ISubsetHandler::assign_subset))
		.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(Vertex*) const>(&ISubsetHandler::get_subset_index))
		.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(Edge*) const>(&ISubsetHandler::get_subset_index))
		.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(Face*) const>(&ISubsetHandler::get_subset_index))
		.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(Volume*) const>(&ISubsetHandler::get_subset_index))

		.add_method("num_subsets", &ISubsetHandler::num_subsets)
		.add_method("get_subset_name", &ISubsetHandler::get_subset_name, "subset name", "subsetIndex")
		.add_method("set_subset_name", &ISubsetHandler::set_subset_name, "", "name#subsetIndex")
		.add_method("get_subset_index", static_cast<int (ISubsetHandler::*)(const char*) const>(
										&ISubsetHandler::get_subset_index), "subsetIndex", "subsetName")
		.add_method("set_default_subset_index", &ISubsetHandler::set_default_subset_index, "", "subsetIndex")
		.add_method("get_default_subset_index", &ISubsetHandler::set_default_subset_index, "subsetIndex", "");
	
//	SubsetHandler
	reg.template add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler", grp)
		.add_constructor()
		.add_method("assign_grid", static_cast<void (SubsetHandler::*)(Grid&)>(&SubsetHandler::assign_grid), "", "g")
		.set_construct_as_smart_pointer(true);


//	MGSubsetHandler
	reg.template add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler", grp)
		.add_constructor()
		.add_method("assign_grid", &MGSubsetHandler::assign_grid, "", "mg")
		.set_construct_as_smart_pointer(true);

//	SurfaceView
	reg.template add_class_<SurfaceView>("SurfaceView", grp)
		.add_method("subset_handler", static_cast<ConstSmartPtr<MGSubsetHandler> (SurfaceView::*)() const>(
										&SurfaceView::subset_handler));

//	subset util
	reg.add_function("AdjustSubsetsForSimulation",
					static_cast<void (*)(SubsetHandler&, bool)>(
					&AdjustSubsetsForSimulation<SubsetHandler>), grp)
		.add_function("AdjustSubsetsForSimulation",
					static_cast<void (*)(MGSubsetHandler&, bool)>(
					&AdjustSubsetsForSimulation<MGSubsetHandler>), grp)
		.add_function("AssignSubsetsByElementType",
					static_cast<void (*)(ISubsetHandler&)>(&AssignSubsetsByElementType))
		.add_function("AssignSubsetColors", &AssignSubsetColors)
		.add_function("AssignSubsetsByLevel", static_cast<void (*)(SubsetHandler&, MultiGrid&)>(&AssignSubsetsByLevel));
}


}//	end of namespace bridge

UG_REGISTRY_DEFINE(RegisterGridBridge_SubsetHandler);

}//	end of namespace ug

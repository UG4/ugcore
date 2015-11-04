// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "grid_bridges.h"
#include "lib_grid/algorithms/subset_util.h"
#include "lib_grid/subset_handler.h"
#include "lib_grid/tools/surface_view.h"
using namespace std;

namespace ug{
namespace bridge{

void RegisterGridBridge_SubsetHandler(Registry& reg, string parentGroup)
{
	string grp = parentGroup;
	
//  ISubsetHandler
	reg.add_class_<ISubsetHandler>("ISubsetHandler", grp)
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
	reg.add_class_<SubsetHandler, ISubsetHandler>("SubsetHandler", grp)
		.add_constructor()
		.add_method("assign_grid", static_cast<void (SubsetHandler::*)(Grid&)>(&SubsetHandler::assign_grid), "", "g")
		.set_construct_as_smart_pointer(true);


//	MGSubsetHandler
	reg.add_class_<MGSubsetHandler, ISubsetHandler>("MGSubsetHandler", grp)
		.add_constructor()
		.add_method("assign_grid", &MGSubsetHandler::assign_grid, "", "mg")
		.set_construct_as_smart_pointer(true);

//	SurfaceView
	reg.add_class_<SurfaceView>("SurfaceView", grp)
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
		.add_function("AssignSubsetColors", &AssignSubsetColors);
}

}//	end of namespace
}//	end of namespace

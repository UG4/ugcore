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

//#include "grid_bridges.h"
#include "registry/registry.h"
#include "common/profiler/profiler.h"
#include "lib_grid/multi_grid.h"
#include "lib_grid/file_io/file_io.h"
#include "lib_grid/file_io/file_io_ugx.h"

using namespace std;

namespace ug{
namespace bridge{

bool LoadGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return LoadGridFromFile(grid, sh, filename);
}

bool LoadGrid(Grid& grid, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return LoadGridFromFile(grid, filename);
}

bool SaveGrid(Grid& grid, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(grid, filename);
}

bool SaveGrid(Grid& grid, ISubsetHandler& sh, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(grid, sh, filename);
}

bool SaveGrid(Grid& grid, const ISubsetHandler& sh, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(grid, *const_cast<ISubsetHandler*>(&sh), filename);
}

bool SaveGridHierarchy(MultiGrid& mg, const char* filename)
{
	PROFILE_FUNC_GROUP("grid");
	return SaveGridToFile(mg, mg.get_hierarchy_handler(), filename);
}


void RegisterGridBridge_FileIO(Registry& reg, string parentGroup)
{
	string grp = parentGroup;

//	UGXFileInfo
	reg.add_class_<UGXFileInfo>("UGXFileInfo", grp)
		.add_constructor()
		.add_method("parse_file", &UGXFileInfo::parse_file, "", "filename")
		.add_method("num_grids", &UGXFileInfo::num_grids)
		.add_method("num_subset_handlers", &UGXFileInfo::num_subset_handlers)
		.add_method("num_subsets", &UGXFileInfo::num_subsets)
		.add_method("grid_name", &UGXFileInfo::grid_name, "grid name", "gridInd")
		.add_method("subset_handler_name", &UGXFileInfo::subset_handler_name, "", "gridInd#shInd")
		.add_method("subset_name", &UGXFileInfo::subset_name, "", "gridInd#shInd#subsetInd")
		.add_method("grid_has_volumes", &UGXFileInfo::grid_has_volumes, "", "gridInd")
		.add_method("grid_has_faces", &UGXFileInfo::grid_has_faces, "", "gridInd")
		.add_method("grid_has_edges", &UGXFileInfo::grid_has_edges, "", "gridInd")
		.add_method("grid_has_vertices", &UGXFileInfo::grid_has_vertices, "", "gridInd")
		.add_method("physical_grid_dimension", &UGXFileInfo::physical_grid_dimension, "", "gridInd")
		.add_method("topological_grid_dimension", &UGXFileInfo::topological_grid_dimension, "", "gridInd")
		.add_method("grid_world_dimension", &UGXFileInfo::grid_world_dimension, "", "gridInd")
		.set_construct_as_smart_pointer(true);

//  GridObject functions
	reg.add_function("LoadGrid", static_cast<bool (*)(Grid&, ISubsetHandler&, const char*)>(&LoadGrid), grp,
			"", "grid#sh#filename")
		.add_function("LoadGrid", static_cast<bool (*)(Grid&, const char*)>(&LoadGrid), grp,
				"", "grid#filename")
		.add_function("SaveGrid", static_cast<bool (*)(Grid&, const ISubsetHandler&, const char*)>(&SaveGrid), grp,
				"", "grid#sh#filename")
		.add_function("SaveGrid", static_cast<bool (*)(Grid&, ISubsetHandler&, const char*)>(&SaveGrid), grp,
				"", "grid#sh#filename")
		.add_function("SaveGrid", static_cast<bool (*)(Grid&, const char*)>(&SaveGrid), grp,
				"", "grid#filename")
//			.add_function("LoadGridObject", &LoadGridObject, grp,
//					"", "go#filename")
//			.add_function("SaveGridObject", &SaveGridObject, grp,
//					"", "go#filename")
		.add_function("SaveGridHierarchy", &SaveGridHierarchy, grp,
				"", "mg#filename")
		.add_function("SaveGridHierarchyTransformed",
					  static_cast<bool (*)(MultiGrid&, ISubsetHandler&, const char*, number)>(
							  &SaveGridHierarchyTransformed),
					  grp, "", "mg#sh#filename#offset")
		.add_function("SaveGridHierarchyTransformed",
					  static_cast<bool (*)(MultiGrid&, const char*, number)>(
							  &SaveGridHierarchyTransformed),
					  grp, "", "mg#filename#offset")
		.add_function("SaveParallelGridLayout", &SaveParallelGridLayout,
				grp, "", "mg#filename#offset")
		.add_function("SaveSurfaceViewTransformed", &SaveSurfaceViewTransformed)
		.add_function("SaveGridLevelToFile", &SaveGridLevelToFile);
}

}//	end of namespace
}//	end of namespace

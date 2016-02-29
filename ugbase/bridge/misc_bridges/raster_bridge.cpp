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

#include <string>
#include "bridge/bridge.h"
#include "bridge/util.h"
#include "bridge/suffix_tag.h"
#include "common/util/raster.h"

using namespace std;

namespace ug{
namespace bridge{

template <class TValue, int TDIM>
static void RegisterRaster(Registry& reg, string name, string grp)
{
	string suffix = GetDimensionSuffix<TDIM>();
	string tag = GetDimensionTag<TDIM>();

	typedef Raster<TValue, TDIM> T;
	string fullName = name + suffix;

	reg.add_class_<T>(fullName, grp)
		.add_method(
			"load_from_asc",
			&T::load_from_asc,
			"", "filename", "Loads the given file and creates the raster accordingly.")
		.add_method(
			"set_num_nodes",
			static_cast<void (T::*)(int, size_t)>(&T::set_num_nodes),
			"", "dim # numNodes", "set the number of nodes for the given dimension.")
		.add_method(
			"num_nodes",
			static_cast<size_t (T::*)(int) const>(&T::num_nodes),
			"numNodes", "dim", "returns the number of nodes for the given dimension.")
		.add_method(
			"create",
			&T::create,
			"", "", "Creates the raster according to the set number of nodes (use 'set_num_nodes').")
		.add_method(
			"set_min_corner",
			static_cast<void (T::*)(int, number)>(&T::set_min_corner),
			"", "dim # coordinate", "set the coordinate of the minimum corner of the raster for the given dimension.")
		.add_method(
			"min_corner",
			static_cast<number (T::*)(int) const>(&T::min_corner),
			"coordinate", "dim", "returns the coordinate of the minimum corner of the raster for the given dimension.")
		.add_method(
			"set_extension",
			static_cast<void (T::*)(int, number)>(&T::set_extension),
			"", "dim # coordinate", "set the extension of the raster for the given dimension.")
		.add_method(
			"extension",
			static_cast<number (T::*)(int) const>(&T::extension),
			"coordinate", "dim", "returns the extension of the raster for the given dimension.")
		.add_method(
			"select_node",
			static_cast<void (T::*)(int, size_t)>(&T::select_node),
			"", "dim # index", "select a node by specifying the index in each dimension.")
		.add_method(
			"selected_node_value",
			&T::selected_node_value,
			"value", "", "returns the value of the selected node (use 'select_node' to select a node).")
		.add_method(
			"set_selected_node_value",
			&T::set_selected_node_value,
			"", "value", "set the value of the selected node (use 'select_node' to select a node).")
		.add_method(
			"set_cursor",
			static_cast<void (T::*)(int, number)>(&T::set_cursor),
			"", "dim # coordinate", "set the coordinate of the cursor for each dimension.")
		.add_method(
			"interpolate_at_cursor",
			&T::interpolate_at_cursor,
			"value", "", "returns the interpolated value (using the given order) at the cursor (use 'set_cursor' to set the cursor).");

	reg.add_class_to_group(fullName, name, tag);
}

void RegisterBridge_Raster(Registry& reg, string parentGroup)
{
	string grp = parentGroup; grp.append("/Util/Raster");

	RegisterRaster<number, 1>(reg, "NumberRaster", parentGroup);
	RegisterRaster<number, 2>(reg, "NumberRaster", parentGroup);
	RegisterRaster<number, 3>(reg, "NumberRaster", parentGroup);
}

}//	end of namespace
}//	end of namespace

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
#include "bridge/suffix_tag.h"
#include "lib_grid/algorithms/deg_layer_mngr.h"
#include "lib_grid/algorithms/extrusion/expand_layers.h"
#include "lib_grid/algorithms/grid_generation/horizontal_layers_mesher.h"

using namespace std;

namespace ug{
namespace bridge{

/**
 * A template function for registering a degenerated layer manager for a
 * specified dimensionality.
 */
template <int dim, typename TRegistry=Registry>
static void RegisterDegeneratedLayerManager(TRegistry& reg, string grp)
{
	string suffix = GetDimensionSuffix<dim>();
	string tag = GetDimensionTag<dim>();

	typedef DegeneratedLayerManager<dim> T;
	string name = string("DegeneratedLayerManager").append(suffix);
	reg.template add_class_<T>(name, grp)
		.template add_constructor<void (*) (SmartPtr<MultiGridSubsetHandler>)>("MultiGridSubsetHandler")
		.add_method("add", static_cast<void (T::*) (const char*)>(&T::add), "Adds degenerated subsets to the manager", "subset(s)")
		.add_method("remove", static_cast<void (T::*) (const char*)>(&T::remove), "Removes subsets from the manager", "subset(s)")
		.add_method("close", static_cast<void (T::*) ()>(&T::close), "Finalizes the fracture manager", "")
		.add_method("contains", static_cast<bool (T::*) (int)>(&T::contains), "Is subset registered in the manager", "subset(s)")
		.add_method("num_subsets", static_cast<size_t (T::*) ()>(&T::num_subsets), "Number of subsets in the manager", "")
		.add_method("subset", static_cast<int (T::*) (size_t)>(&T::subset), "Subset index for a given index in the manager", "index in the manager")
		.add_method("assign_middle_subset", static_cast<int (T::*) (int, int)>(&T::assign_middle_subset), "Assigns the subset index to the middle manifold", "subset index of the layer#subset index for the middle manifold")
		.add_method("assign_middle_subset", static_cast<int (T::*) (int, const char*)>(&T::assign_middle_subset), "Assigns the subset index to the middle manifold", "subset index of the layer#subset name for the middle manifold")
		.add_method("assign_middle_subset", static_cast<int (T::*) (const char*, const char*)>(&T::assign_middle_subset), "Assigns the subset index to the middle manifold", "subset of the layer#subset name for the middle manifold")
		.add_method("init_refiner", static_cast<void (T::*) (SmartPtr<GlobalFracturedMediaRefiner>,bool)>(&T::init_refiner), "Init. refiner", "refiner#as low dim")
		.set_construct_as_smart_pointer(true);
	reg.add_class_to_group(name, "DegeneratedLayerManager", tag);
}

////////////////////////////////////////////////////////////////////////
///	A helper class for ExpandLayers.
/**	This class should never be publicly available, especially since
 * deriving from std::vector is a bad idea (compare 'Effective C++').
 * However, it is very useful in this situation.
 *
 * The class simply extends std::vector<FractureInfo> by an add_layer method.
 */
class ExpandLayersDesc : public std::vector<FractureInfo>
{
	public:
		ExpandLayersDesc() {}

		void add_layer(int subsetInd, int newSubsetInd, number width)
		{
			push_back(FractureInfo(subsetInd, newSubsetInd, width));
		}
};

template <typename TRegistry=Registry>
void RegisterGridBridge_Layers_(TRegistry& reg, string parentGroup)
{
	string grp = parentGroup;
	
#ifdef UG_DIM_2
	RegisterDegeneratedLayerManager<2, TRegistry> (reg, grp);
#endif
#ifdef UG_DIM_3
	RegisterDegeneratedLayerManager<3, TRegistry> (reg, grp);
#endif

	typedef vector<FractureInfo> FracInfoVec;
	reg.template add_class_<FracInfoVec>("FractureInfoVec", grp);

	reg.template add_class_<ExpandLayersDesc, FracInfoVec>("ExpandLayersDesc", grp)
		.add_constructor()
		.add_method("add_layer", &ExpandLayersDesc::add_layer)
		.set_construct_as_smart_pointer(true);

	// \todo: this is uncommented, since in conflict with new vector
	//	handling of registry. Should be adapted.
//		reg.add_function("ExpandLayers2d", &ExpandFractures2d, grp)
//			.add_function("ExpandLayers3d", &ExpandFractures3d, grp);

//	prism-meshing
	reg.template add_class_<RasterLayerDesc>("RasterLayerDesc", grp,
			"Layer Desc for RasterLayers class")
		.template add_constructor<void (*)(const std::string&, number)>(
			"filename#minLayerHeight")
		.add_method("filename", &RasterLayerDesc::filename,
			"filename", "", "Returns the filename of the given layer-desc")
		.add_method("min_height", &RasterLayerDesc::min_height,
			"minHeight", "", "Returns the minimal height of the given layer-desc")
		.set_construct_as_smart_pointer(true);
	
	reg.template add_class_<Heightfield>("Heightfield", grp, "2d raster with number values")
		.add_constructor()
		.add_method("interpolate",
					static_cast<number (Heightfield::*)(number, number) const>(&Heightfield::interpolate),
					"height", "x#y", "returns the height at the given coordinate using piecewise constant interpolation")
		.add_method("interpolate",
					static_cast<number (Heightfield::*)(number, number, int) const>(&Heightfield::interpolate),
					"height", "x#y#order", "returns the height at the given coordinate using the specified interpolation order")
		.add_method("no_data_value", &Heightfield::no_data_value,
					"noDataValue", "", "returns the value which represents an invalid height")
		.add_method("blur", &Heightfield::blur, "", "alpha, iterations",
					"Smoothens the field by adjusting the value of each pixel "
					"towards the average of its neighbours")
		.add_method("eliminate_invalid_cells", &Heightfield::eliminate_invalid_cells,
					"success", "", "eliminates invalid cells by repeatedly filling "
					"those cells with averages of neighboring cells");

	reg.add_function("LoadHeightfieldFromASC", LoadHeightfieldFromASC, grp,
					 "", "heightfield # filename", "Loads a heightfield from the specified file");

	reg.template add_class_<RasterLayers>("RasterLayers", grp, "Stack of 2d raster data.")
		.add_constructor()
		.add_method("load_from_files",
			static_cast<void (RasterLayers::*)(const std::vector<std::string>&, number)>(
				&RasterLayers::load_from_files),
			"", "filenames", "Loads raster data from the specified .asc files. "
			"Specify the bottom layer first and the surface layer last.")
		.add_method("load_from_files",
			static_cast<void (RasterLayers::*)(
				const std::vector<SPRasterLayerDesc>&)>(
					&RasterLayers::load_from_files),
			"", "filenames", "Loads raster data from the specified .asc files. "
			"Specify the bottom layer first and the surface layer last.")
		.add_method("invalidate_flat_cells", &RasterLayers::invalidate_flat_cells, "",
			"min height", "Marks all cells as invalid that belong to a "
			"small lense regarding its horizontal area.")
		.add_method("invalidate_small_lenses", &RasterLayers::invalidate_small_lenses, "",
			"min area", "Marks all cells as invalid which are closer to the "
			"next higher valid cell than the given min height.")
		.add_method("remove_small_holes", &RasterLayers::remove_small_holes, "",
			"max area, min height", "removes small holes by expanding the layer "
			"in those regions to the specified height")
		.add_method("snap_cells_to_higher_layers", &RasterLayers::snap_cells_to_higher_layers, "",
			"min height", "sets invalid or flat cells to the value of the corresponding "
			"cell in the level above.")
		.add_method("eliminate_invalid_cells", &RasterLayers::eliminate_invalid_cells, "success",
			"eliminates invalid cells by filling those cells with averages of neighboring valid cells.")
		.add_method("blur_layers", &RasterLayers::blur_layers, "",
			"alpha # num iterations", "Blurs the values in each layer by averaging between "
			"neighbored cells on the same layer.")
		.add_method("construct_relative_to_global_height_table",
					&RasterLayers::construct_relative_to_global_height_table,
					"", "iterations # alpha",
					"Prepares a table for improved height value reconstructions in no-data-regions.")
		.set_construct_as_smart_pointer(true);

}

}//	end of namespace bridge

UG_REGISTRY_DEFINE(RegisterGridBridge_Layers);

}//	end of namespace ug

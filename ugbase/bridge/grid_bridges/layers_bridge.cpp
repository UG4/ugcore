// created by Sebastian Reiter
// s.b.reiter@gmail.com

#include "grid_bridges.h"
#include "lib_grid/algorithms/extrusion/expand_layers.h"
#include "lib_grid/algorithms/grid_generation/horizontal_layers_mesher.h"

using namespace std;

namespace ug{
namespace bridge{

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


void RegisterGridBridge_Layers(Registry& reg, string parentGroup)
{
	string grp = parentGroup;

	typedef vector<FractureInfo> FracInfoVec;
	reg.add_class_<FracInfoVec>("FractureInfoVec", grp);

	reg.add_class_<ExpandLayersDesc, FracInfoVec>("ExpandLayersDesc", grp)
		.add_constructor()
		.add_method("add_layer", &ExpandLayersDesc::add_layer)
		.set_construct_as_smart_pointer(true);

	// \todo: this is uncommented, since in conflict with new vector
	//	handling of registry. Should be adapted.
//		reg.add_function("ExpandLayers2d", &ExpandFractures2d, grp)
//			.add_function("ExpandLayers3d", &ExpandFractures3d, grp);

//	prism-meshing
	reg.add_class_<RasterLayerDesc>("RasterLayerDesc", grp,
			"Layer Desc for RasterLayers class")
		.add_constructor<void (RasterLayerDesc::*)(const std::string&, number)>(
			"filename#minLayerHeight")
		.add_method("filename", &RasterLayerDesc::filename,
			"filename", "", "Returns the filename of the given layer-desc")
		.add_method("min_height", &RasterLayerDesc::min_height,
			"minHeight", "", "Returns the minimal height of the given layer-desc")
		.set_construct_as_smart_pointer(true);

	reg.add_class_<RasterLayers>("RasterLayers", grp, "Stack of 2d raster data.")
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
		.set_construct_as_smart_pointer(true);

	reg.add_class_<Heightfield>("Heightfield", grp, "2d raster with number values")
		.add_constructor()
		.add_method("interpolate",
					static_cast<number (Heightfield::*)(number, number) const>(&Heightfield::interpolate),
					"height", "x#y", "returns the height at the given coordinate using piecewise constant interpolation")
		.add_method("interpolate",
					static_cast<number (Heightfield::*)(number, number, int) const>(&Heightfield::interpolate),
					"height", "x#y#order", "returns the height at the given coordinate using the specified interpolation order")
		.add_method("no_data_value", &Heightfield::no_data_value,
					"noDataValue", "", "returns the value which represents an invalid height");

	reg.add_function("LoadHeightfieldFromASC", LoadHeightfieldFromASC, grp,
					 "", "heightfield # filename", "Loads a heightfield from the specified file");
}

}//	end of namespace
}//	end of namespace

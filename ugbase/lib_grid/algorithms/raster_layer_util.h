#ifndef __H__UG_raster_layer_util
#define __H__UG_raster_layer_util

#include <utility>
#include <string>
#include <vector>
#include "lib_grid/algorithms/heightfield_util.h"

namespace ug{

struct RasterLayerDesc{
	RasterLayerDesc(const std::string& filename, number minHeight) : 
		m_filename(filename), m_minHeight(minHeight) 	{}

	const std::string& filename() const 	{return m_filename;}
	number min_height() const				{return m_minHeight;}

	private:
		std::string m_filename;
		number		m_minHeight;
};

typedef SmartPtr<RasterLayerDesc>	SPRasterLayerDesc;


class RasterLayers{
	public:
		struct layer_t{
			layer_t() : heightfield(), minHeight(SMALL) {}

			Heightfield heightfield;
			number		minHeight;
		};


		typedef RasterLayerDesc		LayerDesc;
		typedef SPRasterLayerDesc	SPLayerDesc;


	///	loads raster data from a list of .asc files.
	/**	layerDescs.front() represents the bottom of the lowest layer.
	 * layerDescs.top() represents the terrain surface. All data inbetween
	 * is interpreted as sorted layer-bottoms from bottom to top.
	 * \{ */
		void load_from_files(const std::vector<LayerDesc>& layerDescs);

		void load_from_files(const std::vector<SPLayerDesc>& layerDescs);
	/** \} */
	
	///	loads raster data from a list of .asc files.
	/**	filenames.front() represents the bottom of the lowest layer.
	 * filenames.top() represents the terrain surface. All data inbetween
	 * is interpreted as sorted layer-bottoms from bottom to top.
	 *
	 *	\param minLayerHeight	If the height of the layer at a given point is
	 * 							smaller than this value, then the layer is considered
	 *							to be non-existant at this point (i.e. has a hole)*/
		void load_from_files(const std::vector<std::string>& filenames,
							 number minLayerHeight);


		void resize(size_t newSize);
		
		size_t size() const			{return m_layers.size();}
		size_t num_layers() const	{return m_layers.size();}
		bool empty() const			{return m_layers.empty();}

		layer_t& operator[] (size_t i)				{return *m_layers[i];}
		const layer_t& operator[] (size_t i) const	{return *m_layers[i];}

		layer_t& 		layer (size_t i)			{return *m_layers[i];}
		const layer_t&	layer (size_t i) const		{return *m_layers[i];}

		const layer_t&	top () const				{return *m_layers.back();}

		Heightfield&		heightfield (size_t i)			{return layer(i).heightfield;}
		const Heightfield&	heightfield (size_t i) const	{return layer(i).heightfield;}

		void	set_min_height (size_t i, number h)		{layer(i).minHeight = h;}
		number	min_height (size_t i) const				{return layer(i).minHeight;}


	///	invalidates cells in lower levels which are too close to valid cells in higher levels
		void invalidate_flat_cells();

	///	invalidates cells that belong to a small lense regarding its horizontal area
		void invalidate_small_lenses(number minArea);

	///	removes small holes by expanding the layer in those regions to the specified height
		void remove_small_holes(number maxArea);

	///	sets invalid or flat cells to the value of the corresponding cell in the level above
	/** This method is somehow antithetical to 'invalidate_flat_cells', since it reassigns
	 * values to invalid cells which are shadowed by valid cells.*/
		void snap_cells_to_higher_layers();

	///	eliminates invalid cells by filling those cells with averages of neighboring valid cells
		void eliminate_invalid_cells();

	///	smoothens the values in each layer by averaging with neighboured values
		void blur_layers(number alpha, size_t numIterations);

	///	finds the first valid value at the given x-y-coordinate starting at the specified layer moving downwards.
	/** returns a pair containing the layer-index in which the valid value was found (first)
	 * and the value at the given coordinate in that layer. Returns -1 if no such layer
	 * was found.*/
		std::pair<int, number> trace_line_down(const vector2& c, size_t firstLayer) const;

	///	finds the first valid value at the given x-y-coordinate starting at the specified layer moving downwards.
	/** returns a pair containing the layer-index in which the valid value was found (first)
	 * and the value at the given coordinate in that layer. Returns -1 if no such layer
	 * was found.*/
		std::pair<int, number> trace_line_up(const vector2& c, size_t firstLayer) const;

	private:
		std::vector<SmartPtr<layer_t> >	m_layers;
};

}//	end of namespace

#endif	//__H__UG_raster_layer_util

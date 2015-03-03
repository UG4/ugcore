// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_horizontal_layers_mesher
#define __H__UG_horizontal_layers_mesher

#include <utility>
#include <string>
#include <vector>
#include "lib_grid/algorithms/field_util.h"

namespace ug{
	
class RasterLayers{
	public:
		typedef Heightfield layer_t;

	///	loads raster data from a list of .asc files.
	/**	filenames.front() represents the bottom of the lowest layer.
	 * filenames.top() represents the terrain surface. All data inbetween
	 * is interpreted as sorted layer-bottoms from bottom to top.*/
		void load_from_files(const std::vector<std::string>& filenames);
		void resize(size_t newSize);
		
		size_t size() const			{return m_layers.size();}
		size_t num_layers() const	{return m_layers.size();}
		bool empty() const			{return m_layers.empty();}

		layer_t& operator[] (size_t i)				{return *m_layers[i];}
		const layer_t& operator[] (size_t i) const	{return *m_layers[i];}

		layer_t& layer(size_t i)					{return *m_layers[i];}
		const layer_t& layer(size_t i) const		{return *m_layers[i];}

		const layer_t& top() const		{return *m_layers.back();}

	///	invalidates cells in lower levels which are too close to valid cells in higher levels
		void invalidate_flat_cells(number minHeight);

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

void MeshLayerBoundaries(Grid& grid, const RasterLayers& layers,
						 Grid::VertexAttachmentAccessor<AVector3> aaPos,
						 ISubsetHandler* pSH = NULL);

void MeshLayers(Grid& grid, const RasterLayers& layers,
				Grid::VertexAttachmentAccessor<AVector3> aaPos,
				ISubsetHandler* pSH = NULL);

/**	grid has to contain a triangluation of the surface grid of raster-layers.
 * Only x- and y- coordinates of the vertices of the reference triangulation are
 * considered, since all vertices are projected to their respective layers.
 */
void ExtrudeLayers(Grid& grid, const RasterLayers& layers,
				   Grid::VertexAttachmentAccessor<AVector3> aaPos,
				   ISubsetHandler& sh);

}//	end of namespace

#endif	//__H__UG_horizontal_layers_mesher

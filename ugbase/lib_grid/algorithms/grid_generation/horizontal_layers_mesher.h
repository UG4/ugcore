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
		
		size_t size() const	{return m_layers.size();}

		layer_t& operator[] (size_t i)				{return *m_layers[i];}
		const layer_t& operator[] (size_t i) const	{return *m_layers[i];}

		layer_t& layer(size_t i)					{return *m_layers[i];}
		const layer_t& layer(size_t i) const		{return *m_layers[i];}

	///	invalidates cells in lower levels which are too close to valid cells in higher levels
		void invalidate_flat_cells(number minHeight);

	private:
		std::vector<SmartPtr<layer_t> >	m_layers;
};

void MeshLayerBoundaries(Grid& grid, const RasterLayers& layers,
						 Grid::VertexAttachmentAccessor<AVector3> aaPos,
						 ISubsetHandler* pSH = NULL);

void MeshLayers(Grid& grid, const RasterLayers& layers,
				Grid::VertexAttachmentAccessor<AVector3> aaPos,
				ISubsetHandler* pSH = NULL);

}//	end of namespace

#endif	//__H__UG_horizontal_layers_mesher

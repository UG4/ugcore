// created by Sebastian Reiter
// s.b.reiter@gmail.com

#ifndef __H__UG_horizontal_layers_mesher
#define __H__UG_horizontal_layers_mesher

#include "lib_grid/algorithms/raster_layer_util.h"

namespace ug{


void MeshLayerBoundaries(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler* pSH = NULL);

void MeshLayers(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler* pSH = NULL);

/**	grid has to contain a triangluation of the surface grid of raster-layers.
 * Only x- and y- coordinates of the vertices of the reference triangulation are
 * considered, since all vertices are projected to their respective layers.
 * By setting 'allowForTetsAndPyras = true', one will receive less elements. By
 * setting 'allowForTetsAndPyras = false', the resulting mesh will consist of
 * (possibly rather flat) prisms only.
 */
void ExtrudeLayers(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler& sh,
		bool allowForTetsAndPyras);

}//	end of namespace

#endif	//__H__UG_horizontal_layers_mesher

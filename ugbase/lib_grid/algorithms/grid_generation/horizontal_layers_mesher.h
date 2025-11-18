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

#ifndef __H__UG_horizontal_layers_mesher
#define __H__UG_horizontal_layers_mesher

#include "lib_grid/algorithms/raster_layer_util.h"

namespace ug{


void MeshLayerBoundaries(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler* pSH = nullptr);

void MeshLayers(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler* pSH = nullptr);

/**	grid has to contain a triangluation of the surface grid of raster-layers.
 * Only x- and y- coordinates of the vertices of the reference triangulation are
 * considered, since all vertices are projected to their respective layers.
 * By setting 'allowForTetsAndPyras = true', one will receive less elements. By
 * setting 'allowForTetsAndPyras = false', the resulting mesh will consist of
 * (possibly rather flat) prisms only.
 *
 * \param grid		The grid whose surface triangulation will be extended to
 *					a volume mesh.
 * \param layers	A stack of rasters, defining the layered domain.
 * \param aaPos		vertex positions
 * \param sh		Element subsets will be assigned depending on the index of the
 *					layer in which an element is located.
 * \param allowForTetsAndPyras
 *					'true': the resulting mesh will only contain necessary
 *							elements to represent the individual layers.
 *					'false': For each layer in 'layers' a complete layer of
 *							prisms is generated. Individual prisms are then
 *							assigned to subsets depending on which layer contains
 *							their longest edge.
 * \param aRelZOut	(optional) Vertex attachment to which the relative height of
 *					each vertex will be written. This value will equal the layer
 *					index for which a specific vertex has been created.
 */
void ExtrudeLayers(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler& sh,
		bool allowForTetsAndPyras,
		const ANumber* aRelZOut = nullptr);


/**	grid has to contain a triangluation of the surface grid of raster-layers.
 * Only x- and y- coordinates of the vertices of the reference triangulation are
 * considered, since all vertices are projected to their respective layers.
 * The resulting mesh will contain prisms, pyramids, and tetrahedra
 *
 * \param grid		The grid whose surface triangulation will be extended to
 *					a volume mesh.
 * \param layers	A stack of rasters, defining the layered domain.
 * \param aaPos		vertex positions
 * \param sh		Element subsets will be assigned depending on the index of the
 *					layer in which an element is located.
 * \param aRelZOut	(optional) Vertex attachment to which the relative height of
 *					each vertex will be written. This value will equal the layer
 *					index for which a specific vertex has been created.
 */
void ExtrudeLayersMixed (
		Grid& grid, 
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos,
		ISubsetHandler& sh,
		const ANumber* aRelZOut);


///	projects the given (surface-) grid to the specified raster
void ProjectToLayer(
		Grid& grid,
		const RasterLayers& layers,
		int layerIndex,
		Grid::VertexAttachmentAccessor<AVector3> aaPos);

///	horizontally snaps vertices of the given (surface-) grid to the specified raster
void SnapToHorizontalRaster(
		Grid& grid,
		const RasterLayers& layers,
		Grid::VertexAttachmentAccessor<AVector3> aaPos);

}//	end of namespace

#endif
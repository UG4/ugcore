/*
 * expand_layers_arte.h
 *
 *  Created on: 11.07.2024
 *      Author: mknodel
 */

#ifndef UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_EXPAND_LAYERS_ARTE_H_
#define UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_EXPAND_LAYERS_ARTE_H_

namespace ug
{

/**
 * 2 dimensional fracture expansion for finite extensions, using the Arte algorithm
 *
 */
bool ExpandFractures2dArte(Grid& grid, SubsetHandler& sh,
						const std::vector<FractureInfo>& fracInfos,
						bool expandInnerFracBnds, bool expandOuterFracBnds);



}

#endif /* UGBASE_LIB_GRID_ALGORITHMS_EXTRUSION_EXPAND_LAYERS_ARTE_H_ */

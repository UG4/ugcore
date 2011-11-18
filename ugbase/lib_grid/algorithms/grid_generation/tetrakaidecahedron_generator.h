/*
 * tetrakaidekaeder_generator.h
 *
 *  Created on: 18.11.2011
 *      Author: marscher
 */

#ifndef TETRAKAIDEKAEDER_GENERATOR_H_
#define TETRAKAIDEKAEDER_GENERATOR_H_

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug{

void GenerateTetrakaidecahedron(Grid& grid, number& height, number& baseEdgeLength, number& diameter);


/**
 * indsOut: numInds1, ind1_1, ind1_2, ..., numInds2, ind2_1, ind2_2, ...
 *
 * numInds == 4: tetrahedron
 * numInds == 5: pyramid
 * numInds == 6: prism
 * numInds == 8: hexahedron
 */
void GenerateTetrakaidecahedron(std::vector<vector3>& posOut, std::vector<int>& indsOut,
							  number& height, number& baseEdgeLength, number& diameter);

}//	end of namespace

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */

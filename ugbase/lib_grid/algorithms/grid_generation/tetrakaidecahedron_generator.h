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

namespace tkdGenerator{

typedef std::vector<vector3> CoordsArray;
typedef std::vector<int> IndexArray;
typedef const vector3& vec3Ref;

/**
 * \param grid
 * \param height
 * \param baseEdgeLength
 * \param diameter
 */
void GenerateTetrakaidecahedron(Grid& grid, number& height,
		number& baseEdgeLength, number& diameter);

/**
 * indsOut: numInds1, ind1_1, ind1_2, ..., numInds2, ind2_1, ind2_2, ...
 *
 * numInds == 4: tetrahedron
 * numInds == 5: pyramid
 * numInds == 6: prism
 * numInds == 8: hexahedron
 */
void GenerateTetrakaidecahedron(CoordsArray&, IndexArray&,
							  number& height, number& baseEdgeLength, number& diameter);

void createPrism(vec3Ref v1, vec3Ref v2, vec3Ref v3,
				 vec3Ref v4, vec3Ref v5, vec3Ref v6,
				 CoordsArray& posOut, IndexArray& indsOut);

void createTetrahedron(vec3Ref v1, vec3Ref v2, vec3Ref v3, vec3Ref v4,
		CoordsArray& posOut, IndexArray& indsOut);

} // end of namespace tkdGenerator
}//	end of namespace ug

#endif /* TETRAKAIDEKAEDER_GENERATOR_H_ */

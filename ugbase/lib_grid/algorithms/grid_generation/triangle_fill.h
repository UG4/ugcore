//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m05 d14

#ifndef __H__LIB_GRID__TRIANGLE_FILL__
#define __H__LIB_GRID__TRIANGLE_FILL__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/// Fills a 2d-region, bounded by the given poly-chain, with triangles.
bool TriangleFill(std::vector<int>& vTriIndsOut, vector2* polyChain,
					size_t polyChainSize, bool bTriangulateInside = true);

/// Fills a region bounded by the given poly-chain with triangles.
/**
 * This algorithm uses Grid::mark.
 *
 */			
bool TriangleFill(Grid& grid, EdgeBaseIterator edgesBegin,
				EdgeBaseIterator edgesEnd, bool bTriangulateInside = true);


}//	end of namespace

#endif

#ifndef __H__LIB_GRID__TRIANGLE_FILL__
#define __H__LIB_GRID__TRIANGLE_FILL__

#include <vector>
#include "lib_grid/lg_base.h"

namespace ug
{

/// \addtogroup lib_grid_algorithms_grid_generation
///	@{

/// Fills a 2d-region, bounded by the given poly-chain, with triangles.
bool TriangleFill(std::vector<int>& vTriIndsOut, vector2* polyChain,
					size_t polyChainSize, bool bTriangulateInside = true);

/// Fills a region bounded by the given poly-chain with triangles.
/**
 * This algorithm uses Grid::mark.
 *
 */			
bool TriangleFill(Grid& grid, EdgeIterator edgesBegin,
				EdgeIterator edgesEnd, bool bTriangulateInside = true);

/**@}*/ // end of doxygen defgroup command

}//	end of namespace

#endif

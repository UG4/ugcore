//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#ifndef __H__LIB_GRID__EDGE_LENGTH_ADJUSTMENT_EXTENDED__
#define __H__LIB_GRID__EDGE_LENGTH_ADJUSTMENT_EXTENDED__

#include "lib_grid/lg_base.h"
#include "common/ug_config.h"
#include "edge_length_adjustment.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_remeshing
///	@{

struct AdjustEdgeLengthDesc{
	AdjustEdgeLengthDesc();

	number minEdgeLen;
	number maxEdgeLen;
	number approximation;
	number triQuality;
	// number adaptivity;
	bool projectPoints;
};
/**
 * Splits all edges that are too long and collapses edges that
 * are too short.
 *
 * marks: 0 = normal, 1 = crease, 2 = fixed
 */
UG_API
bool AdjustEdgeLength(Grid& grid, SubsetHandler& shMarks,
					  const AdjustEdgeLengthDesc& desc, int numIterations);

/// @}	// end of add_to_group command

}//	end of namespace

#endif

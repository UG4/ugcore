//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#ifndef __H__LIB_GRID__EDGE_LENGTH_ADJUSTMENT__
#define __H__LIB_GRID__EDGE_LENGTH_ADJUSTMENT__

#include "lib_grid/lg_base.h"

namespace ug
{

///	\addtogroup lib_grid_algorithms_remeshing
///	@{

enum RemeshingMarks
{
	RM_NONE = -1,
	RM_CREASE = 0,
	RM_FIXED = 1
};

/**
 * Splits all edges that are too long and collapses edges that
 * are too short.
 *
 * marks: 0 = normal, 1 = crease, 2 = fixed
 */
bool AdjustEdgeLength(Grid& gridOut, SubsetHandler& shOut, SubsetHandler& shMarksOut,
					  Grid& gridIn, SubsetHandler& shIn, SubsetHandler& shMarksIn,
					  number minEdgeLen, number maxEdgeLen, int numIterations,
					  bool projectPoints = true, bool adaptive = true);

/// @}	// end of add_to_group command

}//	end of namespace

#endif

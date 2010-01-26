//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y10 m01 d22

#ifndef __H__LIB_GRID__EDGE_LENGTH_ADJUSTMENT__
#define __H__LIB_GRID__EDGE_LENGTH_ADJUSTMENT__

#include "lib_grid/lg_base.h"

namespace ug
{

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
bool AdjustEdgeLength(Grid& grid, SubsetHandler& shMarks,
					  number minEdgeLen, number maxEdgeLen,
					  int numIterations);


}//	end of namespace

#endif

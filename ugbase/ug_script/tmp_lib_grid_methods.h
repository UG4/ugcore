// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// y10 m01 d25

#ifndef __H__UG__TMP_LIB_GRID_METHODS__
#define __H__UG__TMP_LIB_GRID_METHODS__

////////////////////////////////////////////////////////////////////////
//	You should avoid to include this file in your sources, since it is
//	subject to constant change.
//	Most methods declared here will be moved to lib_grid later on.
////////////////////////////////////////////////////////////////////////

#include "lib_grid/lg_base.h"
#include "lib_grid/algorithms/algorithms.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////
//	SaveSelectedEdgesToObj
bool SaveMarkedEdgesToObj(Grid& grid, const char* filename,
						 SubsetHandler& sh,
						 APosition& aPos = aPosition);
}//	end of namespace

#endif

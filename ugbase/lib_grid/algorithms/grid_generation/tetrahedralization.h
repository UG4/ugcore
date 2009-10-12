//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m09 d17

#ifndef __H__LIB_GRID__TETRAHEDRALIZATION__
#define __H__LIB_GRID__TETRAHEDRALIZATION__

#include "lib_grid/lg_base.h"

namespace ug
{

bool Tetrahedralize(Grid& grid, APosition& aPos = aPosition);

bool Tetrahedralize(Grid& grid, SubsetHandler& sh,
					APosition& aPos = aPosition);

}//	end of namespace

#endif

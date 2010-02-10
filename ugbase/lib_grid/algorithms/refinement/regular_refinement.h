// created by Sebastian Reiter
// s.b.reiter@googlemail.com
//	y10 m02 d09

#ifndef __H__LIB_GRID__REGULAR_REFINER__
#define __H__LIB_GRID__REGULAR_REFINER__

#include "lib_grid/lg_base.h"

namespace ug
{

bool Refine(Grid& grid, Selector& sel, AInt& aInt);

}// end of namespace

#endif

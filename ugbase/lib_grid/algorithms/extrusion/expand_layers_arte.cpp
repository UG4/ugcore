/*
 * expand_layers_arte.cpp
 *
 *  Created on: 11.07.2024
 *      Author: mknodel
 */

#include <boost/function.hpp>

#include "expand_layers.h"
#include "lib_grid/algorithms/geom_obj_util/geom_obj_util.h"
#include "lib_grid/callbacks/callbacks.h"
#include "lib_grid/grid/grid_util.h"

#include <stack>
#include <utility>
#include <vector>
#include <type_traits>
#include <limits>
#include <atomic>
#include <cstddef>
#include <bitset>

#include "support.h"

namespace ug{

bool ExpandFractures2dArte( Grid& grid, SubsetHandler& sh, std::vector<FractureInfo> const & fracInfos,
						    bool expandInnerFracBnds, bool expandOuterFracBnds )
{

	return true;
}
}

#ifndef __H__UG__parallel_refinement__
#define __H__UG__parallel_refinement__

#include "parallel_hanging_node_refiner_multi_grid.h"
#include "parallel_global_fractured_media_refiner.h"
#include "parallel_global_refiner_t.h"

#include "lib_grid/algorithms/refinement/global_multi_grid_refiner.h"
#include "lib_grid/algorithms/refinement/hanging_node_refiner_multi_grid.h"

namespace ug
{

/// \addtogroup lib_grid_parallelization_refinement
/// @{

///	Parallel global refinement for multi-grids
typedef TParallelGlobalRefiner<GlobalMultiGridRefiner>
		ParallelGlobalRefiner_MultiGrid;

///	@}

}//	end of namespace

#endif

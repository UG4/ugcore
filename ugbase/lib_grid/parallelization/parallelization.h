//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y09 m12 d09

#ifndef __H__LIB_GRID__PARALLELIZATION__
#define __H__LIB_GRID__PARALLELIZATION__

#include "parallel_grid_layout.h"
#include "load_balancing.h"
#include "distribution.h"
#include "distributed_grid.h"

#include "copy_policy.h"

#include "parallel_refinement/parallel_refinement.h"
#include "parallelization_util.h"

////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	some doxygen group definitions for various parallelization subgroups follow

////////////////////////////////////////////////////////////////////////
/**
 * \brief parallelization of lib_grid
 *
 * This section contains classes and algorithms that are required to
 * distribute a grid on a parallel architecture, to communicate between
 * parts of a distributed grid and refine a distributed grid.
 *
 * It also contains a sub_group, which contains utility methods for
 * grid distribution and redistribution. 
 *
 * \defgroup lib_grid_parallelization parallelization
 * \ingroup lib_grid
 */

////////////////////////////////////////////////////////////////////////
/**
 * \brief parallel refinement
 *
 * Contains classes and mehtods for parallel refinement.
 *
 * \defgroup lib_grid_parallelization_refinement parallel refinement
 * \ingroup lib_grid_parallelization
 */

////////////////////////////////////////////////////////////////////////
/**
 * \brief parallel distribution
 *
 * Distribution and redistribution of serial and parallel grids.
 *
 * \defgroup lib_grid_parallelization_distribution parallel distribution
 * \ingroup lib_grid_parallelization
 */
 
////////////////////////////////////////////////////////////////////////
/**
 * \brief parallel utilities
 *
 * Contains classes and mehtods that are mostly used internally.
 *
 * \defgroup lib_grid_parallelization_util parallel utilities
 * \ingroup lib_grid_parallelization
 */

#endif

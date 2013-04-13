//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d05

/**
 * \defgroup lib_grid lib_grid
 * \brief The grid library.
 *
 * lib_grid not only defines basic types like ug::Grid and ug::MultiGrid,
 * it also features tools such as the ug::Selector, ug::MGSelector,
 * ug::SubsetHandler and ug::MGSubsetHandler as well as a vast collection
 * of methods that work on those types. Those can mostly be found in the
 * \ref lib_grid_algorithms section.
 *
 */
 
#include "grid/grid.h"
#include "grid/grid_util.h"
#include "geometric_objects/geometric_objects.h"
#include "common_attachments.h"
#include "selector.h"
#include "subset_handler.h"
#include "tools/surface_view.h"
#include "file_io/file_io.h"
#include "algorithms/algorithms.h"
#include "multi_grid.h"

#ifdef UG_PARALLEL
	#include "lib_grid/parallelization/parallelization.h"
#endif

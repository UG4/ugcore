//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m11 d05

#include "grid/grid.h"
#include "grid/grid_util.h"
#include "geometric_objects/geometric_objects.h"
#include "common_attachments.h"
#include "selector.h"
#include "subset_handler.h"
#include "file_io/file_io.h"
#include "algorithms/algorithms.h"
#include "multi_grid.h"

#ifdef UG_PARALLEL
	#include "lib_grid/algorithms/parallelization/parallelization.h"
#endif

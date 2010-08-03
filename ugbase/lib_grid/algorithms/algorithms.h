//	created by Sebastian Reiter
//	s.b.reiter@googlemail.com
//	y08 m12 d04

#ifndef __H__LIB_GRID__ALGORITHMS__
#define __H__LIB_GRID__ALGORITHMS__

#include "trees/kd_tree_static.h"
#include "trees/octree.h"
#include "geom_obj_util/geom_obj_util.h"
#include "grid_generation/grid_generation.h"
#include "refinement/hanging_node_refiner.h"
#include "refinement/regular_refinement.h"
#include "refinement/multi_grid_refiner.h"
#include "refinement/global_multi_grid_refiner.h"
#include "extrusion/extrusion.h"
#include "subdivision/subdivision_loop.h"
#include "extruder_util.h"
#include "subset_util.h"
#include "attachment_util.h"
#include "serialization.h"
#include "selection_util.h"
#include "grid_statistics.h"
#include "multi_grid_util.h"
#include "graph/graph.h"
#include "remeshing/edge_length_adjustment.h"
#include "remeshing/grid_adaption.h"


////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////
//	some doxygen group definitions for various algorithm subgroups follow

////////////////////////////////////////////////////////////////////////
/**
 * \brief contains a variety of useful algorithms.
 *
 * The algorithms section contains a collection of algorithms that can be
 * used to manipulate a grid, obtain information about a grid. It also
 * contains many little helper functions that ease programming a lot.
 *
 * \defgroup lib_grid_algorithms algorithms
 * \ingroup lib_grid
 */

////////////////////////////////////////////////////////////////////////
/**
 * \brief refinement classes and algorithms
 *
 * The refinement section contains classes that allow to refine a
 * ug::Grid or a ug::MultiGrid with different methods.
 * An adaptive refiner that generates a grid that contains hanging vertices
 * (ug::HangingVertex) and constrained and constraining edges
 * (ug::ConstrainedEdge, ug::ConstrainingEdge) is featured as well
 * as a refiner that builds a regular closure.
 *
 * Please note that parallel refinement is adressed in the section
 * \ref lib_grid_parallelization_refinement.
 *
 * \defgroup lib_grid_algorithms_refinement refinement
 * \ingroup lib_grid_algorithms
 */

#endif

/*
 * Copyright (c) 2009-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__LIB_GRID__ALGORITHMS__
#define __H__LIB_GRID__ALGORITHMS__

#include "debug_util.h"
#include "trees/kd_tree_static.h"
#include "trees/octree.h"
#include "geom_obj_util/geom_obj_util.h"
#include "grid_generation/grid_generation.h"
#include "extrusion/extrusion.h"
#include "subdivision/subdivision_loop.h"
#include "extruder_util.h"
#include "subset_util.h"
#include "subset_dim_util.h"
#include "attachment_util.h"
#include "serialization.h"
#include "selection_util.h"
#include "grid_statistics.h"
#include "multi_grid_util.h"
#include "graph/graph.h"
#include "remeshing/edge_length_adjustment.h"
#include "remeshing/grid_adaption.h"
// #include "remeshing/delaunay_triangulation.h"
#include "duplicate.h"
#include "volume_calculation.h"

// hanging_node_refiner_2d_irn.h is currently not maintained.
//#include "refinement/hanging_node_refiner_2d_irn.h"


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

////////////////////////////////////////////////////////////////////////
/**
 * \brief remeshing algorithms
 *
 * Remeshing algorithms change the topology of a grid in order to
 * e.g. adapt it to a given shape or to optimize its triangle aspect ratios.
 *
 * \defgroup lib_grid_algorithms_remeshing remeshing
 * \ingroup lib_grid_algorithms
 */

////////////////////////////////////////////////////////////////////////
/**
 * \brief tree structures
 *
 * Trees allow to perform space-partitioning on a grid. This is important
 * if collision-tests, vertex-projection or clustering shall be performed
 * with maximal performance.
 *
 * \defgroup lib_grid_algorithms_trees trees
 * \ingroup lib_grid_algorithms
 */

////////////////////////////////////////////////////////////////////////
/**
 * \brief subdivision classes and algorithms
 *
 * Subdivision is a special form of refinement, where smooth surfaces
 * are described through repeated refinement and application of position masks.
 *
 * \defgroup lib_grid_algorithms_refinement_subdivision subdivision
 * \ingroup lib_grid_algorithms_refinement
 */

////////////////////////////////////////////////////////////////////////
/**
 * \brief grid generation algorithms
 *
 * Those algorithms allow to construct a triangulation or a tetrahedralization
 * of a given domain.
 *
 * \defgroup lib_grid_algorithms_grid_generation grid generation
 * \ingroup lib_grid_algorithms
 */


#endif

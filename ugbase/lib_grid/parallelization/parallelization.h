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

#ifndef __H__LIB_GRID__PARALLELIZATION__
#define __H__LIB_GRID__PARALLELIZATION__

#include "parallel_grid_layout.h"
#include "load_balancing.h"
#include "distribution.h"
#include "distributed_grid.h"

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

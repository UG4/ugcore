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

////////////////////////////////////////////////////////////////////////
//	This file includes the most basic components of libGrid.
//	You should include this file instead of lib_grid.h whenever
//	you're creating a file that is contained in libGrid.
////////////////////////////////////////////////////////////////////////

#ifndef __H__LIB_GRID__LG_BASE__
#define __H__LIB_GRID__LG_BASE__

#include "common_attachments.h"
#include "lib_grid_messages.h"
#include "multi_grid.h"
#include "selector.h"
#include "subset_handler.h"
#include "grid/grid.h"
#include "grid/grid_util.h"
#include "grid/neighborhood.h"
#include "grid_objects/grid_objects.h"
#include "iterators/lg_for_each.h"
#include "tools/marker_points.h"
#include "tools/bool_marker.h"

#endif

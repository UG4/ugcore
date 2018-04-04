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

#ifndef __H__LIBGRID__FILE_IO_BIN__
#define __H__LIBGRID__FILE_IO_BIN__

#include "lib_grid/grid/grid.h"
#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/selector_interface.h"
#include "lib_grid/common_attachments.h"

namespace ug
{

class ProjectionHandler;

/**
 * Saves a grid to LibGridBinary-format.
 * Awaits a list of subset-handler-pointers and the number
 * of subset-handlers that shall be written.
 */
bool SaveGridToLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ISelector** ppSel, int numSels,
				   ProjectionHandler* pPH = NULL,
				   APosition aPos = aPosition);

bool SaveGridToLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ProjectionHandler* pPH = NULL,
				   APosition aPos = aPosition);


/**
 * Loads a grid from LibGridBinary-format.
 * Awaits a list of subset-handler-pointers and the number
 * of subset-handlers that shall be read.
 * Make sure that all passed subset-handlers are already registered
 * at the grid.
 */
bool LoadGridFromLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ISelector** ppSel, int numSels,
				   ProjectionHandler* pPH = NULL,
				   APosition aPos = aPosition);


bool LoadGridFromLGB(Grid& grid, const char* filename,
				   ISubsetHandler** ppSH, int numSHs,
				   ProjectionHandler* pPH = NULL,
				   APosition aPos = aPosition);

}//	end of namespace

#endif

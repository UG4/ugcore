/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_GRID__SUBSET_DIM_UTIL__
#define __H__UG__LIB_GRID__SUBSET_DIM_UTIL__


#include "lib_grid/tools/subset_handler_interface.h"
#include "lib_grid/tools/subset_handler_multi_grid.h"
#include "lib_grid/tools/subset_handler_grid.h"

namespace ug{

////////////////////////////////////////////////////////////////////////
/// returns if a subset is a regular grid
/**
 * This function returns if a subset contains constrained/constraining elements
 * such as hanging vertices, contrained edges/faces. In this case, the subset
 * does not form a regular grid.
 *
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		true		if subset is regular grid
 * 				false 		if subset is non-regular grid
 */
bool SubsetIsRegularGrid(const SubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
/// returns if a subset is a regular grid
/**
 * This function returns if a subset contains constrained/constraining elements
 * such as hanging vertices, contrained edges/faces. In this case, the subset
 * does not form a regular grid.
 *
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		true		if subset is regular grid
 * 				false 		if subset is non-regular grid
 */
bool SubsetIsRegularGrid(const MGSubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
/// returns if a subset is a regular grid
/**
 * This function returns if a subset contains constrained/constraining elements
 * such as hanging vertices, contrained edges/faces. In this case, the subset
 * does not form a regular grid.
 *
 *
 * \param[in]	sh			SubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		true		if subset is regular grid
 * 				false 		if subset is non-regular grid
 */
bool SubsetIsRegularGrid(const ISubsetHandler& sh, int si);

/// abbreviations for return types
enum {DIM_SUBSET_EMPTY_GRID = -1};

////////////////////////////////////////////////////////////////////////
///	Returns the dimension of geometric objects, that are contained in the subset
/**
 * This function returns the dimension of the subset. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained in the subset
 * If a InterfaceCommunicator is passed, the highest dimension within all
 * procs in the ProcessCommunicator is returned.
 *
 * \param[in]	sh			ISubsetHandler
 * \param[in]	si			Subset Index
 *
 * \return		dimension					Dimension of Subset
 * 				DIM_SUBSET_EMPTY_GRID		if empty Grid given
 */
int DimensionOfSubset(const ISubsetHandler& sh, int si);

////////////////////////////////////////////////////////////////////////
///	Returns the dimension of geometric objects, that are contained in the subset handler
/**
 * This function returns the dimension of the subsets. The dimension is simply
 * defined to be the highest reference dimension of all geometric objects
 * contained the union of all subset
 * If a InterfaceCommunicator is passed, the highest dimension within all
 * procs in the ProcessCommunicator is returned.
 *
 * \param[in]	sh			ISubsetHandler
 *
 * \return		dimension					Dimension of Subset
 * 				DIM_SUBSET_EMPTY_GRID		if empty Grid given
 */
int DimensionOfSubsets(const ISubsetHandler& sh);


} // end namespace ug

#endif /* __H__UG__LIB_GRID__SUBSET_DIM_UTIL__ */

/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__UG__LIB_DISC__COMMON__GROUPS_UTIL__
#define __H__UG__LIB_DISC__COMMON__GROUPS_UTIL__

#include <vector>

#include "lib_disc/common/function_group.h"
#include "lib_disc/dof_manager/function_pattern.h"

namespace ug {

////////////////////////////////////////////////////////////////////////////////
//	Index Mapping
////////////////////////////////////////////////////////////////////////////////

/**
 * Creates a function index mapping that maps all local indices from the one
 * Function Group to the other. Make sure that the first
 * Function Group is contained in the second and both function groups are based
 * on the same function pattern - otherwise an exception is thrown.
 */
void CreateFunctionIndexMapping(FunctionIndexMapping& map,
                                const FunctionGroup& grpFrom,
                                const FunctionGroup& grpTo);

/**
 * Creates a function index mapping that maps all local indices from the one
 * Function Group to the corresponding indices in a function pattern.
 * Make sure that the Function Group is contained function pattern - otherwise
 * an exception is thrown.
 */
void CreateFunctionIndexMapping(FunctionIndexMapping& map,
                                const FunctionGroup& grpFrom,
                                ConstSmartPtr<FunctionPattern> fctPattern);

/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain nullptr)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
void CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const std::vector<const FunctionGroup*>& vFctGrp,
                                 bool sortFct = false);

} // end namespace ug

#endif
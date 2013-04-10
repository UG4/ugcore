/*
 * groups_util.h
 *
 *  Created on: 18.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__COMMON__GROUPS_UTIL__
#define __H__UG__LIB_DISC__COMMON__GROUPS_UTIL__

#include <vector>
#include <string>
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/dof_manager/function_pattern.h"

namespace ug{

////////////////////////////////////////////////////////////////////////////////
//	Subset Group
////////////////////////////////////////////////////////////////////////////////

/**
 * Returns if dimension is the same in all subsets of the subset group
 * @param subsetGroup	subset group that is checked
 * @returns true if dimension is the same in all subsets, else false
 */
bool SameDimensionsInAllSubsets(const SubsetGroup& subsetGroup);

/**
 * Removes all subsets from the subset group that have a lower dimension than the
 * highest dimension contained in the subset group.
 * @param subsetGroup 	subset group that is modified
 */
void RemoveLowerDimSubsets(SubsetGroup& subsetGroup);

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
                                const FunctionPattern& fctPattern);

/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain NULL)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
void CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const std::vector<const FunctionGroup*>& vFctGrp,
                                 bool sortFct = false);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__COMMON__GROUPS_UTIL__ */

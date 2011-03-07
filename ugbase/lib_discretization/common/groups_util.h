/*
 * groups_util.h
 *
 *  Created on: 18.10.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__COMMON__GROUPS_UTIL__
#define __H__LIB_DISCRETIZATION__COMMON__GROUPS_UTIL__

#include <vector>
#include <string>
#include "./subset_group.h"
#include "./function_group.h"
#include "../dof_manager/function_pattern.h"

namespace ug{

// predeclaration
class FunctionPattern;

/**
 * Passing a string of subset names separated by ',' this function returns
 * a subset group containing the subsets.
 */
bool
ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const ISubsetHandler& sh,
                           const char* subsets, const char separator = ',');

/**
 * Passing a string of subset names separated by ',' this function returns
 * a subset group containing the subsets.
 */
bool
ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const FunctionPattern& pattern,
                           const char* subsets, const char separator = ',');

/**
 * Passing a string of function names separated by ',' this function returns
 * a function group containing the subsets.
 */
bool
ConvertStringToFunctionGroup(FunctionGroup&functionGroup, const FunctionPattern& pattern,
                             const char* functions, const char separator = ',');

/**
 * Creates a function index mapping that maps all local indices from the one
 * Function Group to the other. This is of coarse only possible if the first
 * Function Group is contained in the second.
 */
bool
CreateFunctionIndexMapping(FunctionIndexMapping& map,
                           const FunctionGroup& grpFrom,
                           const FunctionGroup& grpTo);


/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain NULL)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
bool CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const std::vector<FunctionGroup*>& vFctGrp,
                                 bool sortFct = false);

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__COMMON__GROUPS_UTIL__ */

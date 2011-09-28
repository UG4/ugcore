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
 * Passing a vector of subset names this function returns
 * a subset group containing the subsets.
 */
bool ConvertStringToSubsetGroup(	SubsetGroup& subsetGroup,
                                	const ISubsetHandler& sh,
									const std::vector<std::string>& vSS);

/**
 * Passing a string of subset names separated by ',' this function returns
 * a subset group containing the subsets.
 */
bool
ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const FunctionPattern& pattern,
                           const char* subsets, const char separator = ',');

/**
 * Passing a string of function names separated by ',' this function returns
 * a function group containing the functions.
 */
bool
ConvertStringToFunctionGroup(FunctionGroup&functionGroup, const FunctionPattern& pattern,
                             const char* functions, const char separator = ',');

/**
 * Passing a vector of function names this function returns
 * a function group containing the functions.
 */
bool ConvertStringToFunctionGroup(	FunctionGroup& functionGroup,
                                  	const FunctionPattern& pattern,
									const std::vector<std::string>& vFct);

/**
 * Creates a function index mapping that maps all local indices from the one
 * Function Group to the other. This is only possible if the first
 * Function Group is contained in the second.
 */
bool
CreateFunctionIndexMapping(FunctionIndexMapping& map,
                           const FunctionGroup& grpFrom,
                           const FunctionGroup& grpTo);

/**
 * Creates a function index mapping that maps all local indices from the one
 * Function Group to the other. The second function group is here contained in
 * the first, thus not all mappings are defined, but only those where the index
 * of the first group is contained in the range of the second.
 */
bool
CreateFunctionIndexMappingInverse(FunctionIndexMapping& map,
                                  const FunctionGroup& grpFromLarge,
                                  const FunctionGroup& grpToSmall);


/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain NULL)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
bool CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const std::vector<const FunctionGroup*>& vFctGrp,
                                 bool sortFct = false);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__COMMON__GROUPS_UTIL__ */

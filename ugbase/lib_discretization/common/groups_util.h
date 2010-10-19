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

bool ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const FunctionPattern& pattern,
								const char* subsets, std::string separator = ",");

bool ConvertStringToFunctionGroup(	FunctionGroup&functionGroup, const FunctionPattern& pattern,
									const char* functions, std::string separator = ",");


} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__COMMON__GROUPS_UTIL__ */

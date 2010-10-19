/*
 * groups_util.cpp
 *
 *  Created on: 18.10.2010
 *      Author: andreasvogel
 */

#include "groups_util.h"
#include "common/string_util.h"

namespace ug{

bool ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const FunctionPattern& pattern,
								const char* subsets, std::string separator)
{
//	get strings
	std::string subsetString = std::string(subsets);

//	set underlying subsethandler Subset Group
	subsetGroup.set_subset_handler(*(pattern.get_subset_handler()));

//	tokenize strings and select subsets
	std::vector<std::string> tokens;
	TokenizeString(subsetString, tokens, separator);

	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);

		if(!subsetGroup.add_subset(tokens[i].c_str()))
		{
			UG_LOG("Name of subset ('" << tokens[i] << "') not found in Subset Handler.\n");
			return false;
		}
	}

//	we're done
	return true;
}


bool ConvertStringToFunctionGroup(	FunctionGroup& functionGroup, const FunctionPattern& pattern,
									const char* functions, std::string separator)
{
//	get strings
	std::string fctString = std::string(functions);

//	create Function Group and Subset Group
	functionGroup.set_function_pattern(pattern);

//	tokenize strings and select functions
	std::vector<std::string> tokens;
	TokenizeString(fctString, tokens, ",");

	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);

		if(!functionGroup.add_function(tokens[i].c_str()))
		{
			UG_LOG("Name of function ('" << tokens[i] << "') not found in Function Pattern.\n");
			return false;
		}
	}

	return true;
}

} // end namespace ug

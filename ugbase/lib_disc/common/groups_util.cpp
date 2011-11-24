/*
 * groups_util.cpp
 *
 *  Created on: 18.10.2010
 *      Author: andreasvogel
 */

#include "groups_util.h"
#include "common/util/string_util.h"

namespace ug{

bool ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const FunctionPattern& pattern,
								const char* subsets, const char separator)
{
//	forward request
	return ConvertStringToSubsetGroup(subsetGroup,
	                                  pattern.subset_handler(),
	                                  subsets, separator);
}

bool ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const ISubsetHandler& sh,
								const char* subsets, const char separator)
{
//	get strings
	std::string subsetString = std::string(subsets);

//	set underlying subsethandler Subset Group
	subsetGroup.set_subset_handler(sh);

//	tokenize strings and select subsets
	std::vector<std::string> tokens;
	TokenizeString(subsetString, tokens, separator);

	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);

		if(!subsetGroup.add(tokens[i].c_str()))
		{
			UG_LOG("Name of subset ('" << tokens[i] << "') not found in Subset Handler.\n");
			return false;
		}
	}

//	we're done
	return true;
}

bool ConvertStringToSubsetGroup(	SubsetGroup& subsetGroup,
                                	const ISubsetHandler& sh,
									const std::vector<std::string>& vSS)
{
//	create Function Group and Subset Group
	subsetGroup.set_subset_handler(sh);

//	tokenize strings and select functions
	for(size_t i = 0; i < vSS.size(); ++i)
		if(!subsetGroup.add(vSS[i].c_str()))
		{
			UG_LOG("ERROR in 'ConvertStringToSubsetGroup': Name of subset ('"
					<< vSS[i] << "') not found in Subset Handler.\n");
			return false;
		}

//	done
	return true;
}


bool ConvertStringToFunctionGroup(	FunctionGroup& functionGroup, const FunctionPattern& pattern,
									const char* functions, const char separator)
{
//	get strings
	std::string fctString = std::string(functions);

//	create Function Group and Subset Group
	functionGroup.set_function_pattern(pattern);

//	tokenize strings and select functions
	std::vector<std::string> tokens;
	TokenizeString(fctString, tokens, ',');

	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);

		if(!functionGroup.add(tokens[i].c_str()))
		{
			UG_LOG("Name of function ('" << tokens[i] << "') not found in"
					" Function Pattern.\n");
			return false;
		}
	}

	return true;
}

bool ConvertStringToFunctionGroup(	FunctionGroup& functionGroup,
                                  	const FunctionPattern& pattern,
									const std::vector<std::string>& vFct)
{
//	create Function Group and Subset Group
	functionGroup.set_function_pattern(pattern);

//	tokenize strings and select functions
	for(size_t i = 0; i < vFct.size(); ++i)
		if(!functionGroup.add(vFct[i].c_str()))
		{
			UG_LOG("ERROR in 'ConvertStringToFunctionGroup': Name of function ('"
					<< vFct[i] << "') not found in Function Pattern.\n");
			return false;
		}

//	done
	return true;
}


bool
CreateFunctionIndexMapping(FunctionIndexMapping& map,
                           const FunctionGroup& grpFromSmall,
                           const FunctionGroup& grpToLarge)
{
//	clear map
	map.clear();

//	check that from group is contained in to group
	if(!grpToLarge.contains(grpFromSmall))
	{
		UG_LOG("ERROR in 'CreateFunctionIndexMapping': smaller FunctionGroup "
				<< grpFromSmall << " is not contained in larger Group " <<
				grpToLarge<<". Cannot create Mapping.\n");
		return false;
	}

//	loop all functions on grpFrom
	for(size_t from = 0; from < grpFromSmall.num_fct(); ++from)
	{
	//	get unique id of function
		const size_t uniqueID = grpFromSmall[from];

	//	find unique id of function in grpTo
		const size_t locIndex = grpToLarge.local_index(uniqueID);

	//	set mapping
		map.add(from, locIndex);
	}

//	we're done
	return true;
}

bool
CreateFunctionIndexMappingInverse(FunctionIndexMapping& map,
                                  const FunctionGroup& grpFromLarge,
                                  const FunctionGroup& grpToSmall)
{
//	clear map
	map.clear();

//	check that from group is contained in to group
	if(!grpFromLarge.contains(grpToSmall)) return false;

//	loop all functions on grpFrom
	for(size_t to = 0; to < grpToSmall.num_fct(); ++to)
	{
	//	get unique id of function
		const size_t uniqueID = grpToSmall[to];

	//	find unique id of function in grpTo
		const size_t locIndex = grpFromLarge.local_index(uniqueID);

	//	set mapping
		map.add(locIndex, to);
	}

//	we're done
	return true;
}


/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain NULL)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
bool CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const std::vector<const FunctionGroup*>& vFctGrp,
                                 bool sortFct)
{
//	clear group
	fctGrp.clear();

//	if empty, nothing to do
	if(vFctGrp.empty()) return true;

//	set underlying subsetHandler
	size_t grp = 0;
	for(; grp < vFctGrp.size(); ++grp)
	{
		if(vFctGrp[grp] == NULL) continue;

		const FunctionPattern* pFctPat = vFctGrp[grp]->get_function_pattern();
		if(pFctPat == NULL)
		{
			UG_LOG("ERROR in 'CreateUnionOfFunctionGroups': Function group "
					<<grp<<" has NULL as underlying FunctionPattern.\n");
			return false;
		}
		fctGrp.set_function_pattern(*pFctPat);
		break;
	}

//	if no function group given
	if(grp == vFctGrp.size()) return true;

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vFctGrp.size(); ++i)
	{
	//	add subset group of elem disc
		if(vFctGrp[i] != NULL)
			if(!fctGrp.add(*vFctGrp[i]))
			{
				UG_LOG("ERROR in 'CreateUnionOfFunctionGroups': Cannot add "
						"functions of the Function Group "<< i << ".\n");
				return false;
			}
	}

//	sort iff required
	if(sortFct) fctGrp.sort();

//	we're done
	return true;
}


} // end namespace ug

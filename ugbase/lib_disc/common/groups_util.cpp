/*
 * groups_util.cpp
 *
 *  Created on: 18.10.2010
 *      Author: andreasvogel
 */

#include "groups_util.h"
#include "common/util/string_util.h"

#include <algorithm>
#include <limits>

namespace ug{

bool SameDimensionsInAllSubsets(const SubsetGroup& subsetGroup)
{
//	compute maximum
	int max = std::numeric_limits<int>::min();
	for(size_t s = 0; s < subsetGroup.size(); ++s)
		max = std::max(subsetGroup.dim(s), max);

//	check
	for(size_t s = 0; s < subsetGroup.size(); ++s)
		if(subsetGroup.dim(s) < max)
			return false;

//	same dimension in all subsets
	return true;
}

void RemoveLowerDimSubsets(SubsetGroup& subsetGroup)
{
//	compute maximum
	int max = std::numeric_limits<int>::min();
	for(size_t s = 0; s < subsetGroup.size(); ++s)
		max = std::max(subsetGroup.dim(s), max);

//	check
	size_t s = 0;
	while(s < subsetGroup.size())
	{
		if(subsetGroup.dim(s) < max)
		{
			// remove and start again
			subsetGroup.remove(subsetGroup[s]);
			s = 0;
		}
		else
		{
			// next
			++s;
		}
	}
}


void ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, const FunctionPattern& pattern,
								const char* subsets, const char separator)
{
//	forward request
	ConvertStringToSubsetGroup(subsetGroup,
	                           pattern.subset_handler(),
	                           subsets, separator);
}


void ConvertStringToSubsetGroup(SubsetGroup& subsetGroup, ConstSmartPtr<ISubsetHandler> pSH,
								const char* subsets, const char separator)
{
//	set underlying subsethandler Subset Group
	subsetGroup.set_subset_handler(pSH);

//	forward
	ConvertStringToSubsetGroup(subsetGroup, subsets, separator);
}

void ConvertStringToSubsetGroup(SubsetGroup& subsetGroup,
								const char* subsets, const char separator)
{
//	get strings
	std::string subsetString = std::string(subsets);

//	tokenize strings and select subsets
	std::vector<std::string> tokens;
	TokenizeString(subsetString, tokens, separator);

	for(size_t i = 0; i < tokens.size(); ++i)
	{
		RemoveWhitespaceFromString(tokens[i]);

		try{
			subsetGroup.add(tokens[i].c_str());
		}UG_CATCH_THROW("Name of subset ('" << tokens[i] <<
		                "') not found in Subset Handler.");
	}
}

void ConvertStringToSubsetGroup(	SubsetGroup& subsetGroup,
                                	ConstSmartPtr<ISubsetHandler> pSH,
									const std::vector<std::string>& vSS)
{
//	create Function Group and Subset Group
	subsetGroup.set_subset_handler(pSH);

//	tokenize strings and select functions
	for(size_t i = 0; i < vSS.size(); ++i){
		try{
			subsetGroup.add(vSS[i].c_str());
		}UG_CATCH_THROW("Name of subset ('" << vSS[i] <<
		                "') not found in Subset Handler.");
	}
}


void ConvertStringToFunctionGroup(	FunctionGroup& functionGroup, const FunctionPattern& pattern,
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

		try{
			functionGroup.add(tokens[i].c_str());
		}UG_CATCH_THROW("Name of function ('" << tokens[i] << "') not found in"
					" Function Pattern.");
	}
}

void ConvertStringToFunctionGroup(	FunctionGroup& functionGroup,
                                  	const FunctionPattern& pattern,
									const std::vector<std::string>& vFct)
{
//	create Function Group and Subset Group
	functionGroup.set_function_pattern(pattern);

//	tokenize strings and select functions
	for(size_t i = 0; i < vFct.size(); ++i)
	{
		try{
			functionGroup.add(vFct[i].c_str());
		}UG_CATCH_THROW("Name of function ('" << vFct[i] << "') not found "
		                "in Function Pattern.");
	}
}


void
CreateFunctionIndexMapping(FunctionIndexMapping& map,
                           const FunctionGroup& grpFromSmall,
                           const FunctionGroup& grpToLarge)
{
//	clear map
	map.clear();

//	check that from group is contained in to group
	if(!grpToLarge.contains(grpFromSmall))
		UG_THROW("CreateFunctionIndexMapping: smaller FunctionGroup "
				<< grpFromSmall << " is not contained in larger Group " <<
				grpToLarge<<". Cannot create Mapping.");

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
}

void
CreateFunctionIndexMappingInverse(FunctionIndexMapping& map,
                                  const FunctionGroup& grpFromLarge,
                                  const FunctionGroup& grpToSmall)
{
//	clear map
	map.clear();

//	check that from group is contained in to group
	if(!grpFromLarge.contains(grpToSmall))
		UG_THROW("CreateFunctionIndexMapping: smaller FunctionGroup "
				<< grpToSmall << " is not contained in larger Group " <<
				grpFromLarge<<". Cannot create Mapping.");

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
}


/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain NULL)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
void CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const std::vector<const FunctionGroup*>& vFctGrp,
                                 bool sortFct)
{
//	clear group
	fctGrp.clear();

//	if empty, nothing to do
	if(vFctGrp.empty()) return;

//	set underlying subsetHandler
	size_t grp = 0;
	for(; grp < vFctGrp.size(); ++grp)
	{
		if(vFctGrp[grp] == NULL) continue;

		const FunctionPattern* pFctPat = vFctGrp[grp]->function_pattern();
		if(pFctPat == NULL)
			UG_THROW("CreateUnionOfFunctionGroups: Function group "
					<<grp<<" has NULL as underlying FunctionPattern.");

		fctGrp.set_function_pattern(*pFctPat);
		break;
	}

//	if no function group given
	if(grp == vFctGrp.size()) return;

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vFctGrp.size(); ++i)
	{
	//	add subset group of elem disc
		if(vFctGrp[i] != NULL)
		{
			try{
				fctGrp.add(*vFctGrp[i]);
			}UG_CATCH_THROW("Cannot add functions of the Function Group "<< i << ".");
		}
	}

//	sort iff required
	if(sortFct) fctGrp.sort();
}


} // end namespace ug

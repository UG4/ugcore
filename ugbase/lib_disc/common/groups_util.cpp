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

#include "groups_util.h"
#include "common/util/string_util.h"

#include <algorithm>
#include <limits>

using namespace std;

namespace ug{

void
CreateFunctionIndexMapping(FunctionIndexMapping& map,
                           const FunctionGroup& grpFromSmall,
                           const FunctionGroup& grpToLarge)
{
//	clear map
	map.clear();

//	check that groups are based on same function pattern
	if(grpToLarge.function_pattern() != grpFromSmall.function_pattern())
		UG_THROW("CreateFunctionIndexMapping: groups are not based on same "
				"function pattern, thus no mapping can be created.");

//	check that "from group" is contained in "to group"
	if(!grpToLarge.contains(grpFromSmall))
		UG_THROW("CreateFunctionIndexMapping: smaller FunctionGroup "
				<< grpFromSmall << " is not contained in larger Group " <<
				grpToLarge<<". Cannot create Mapping.");

//	loop all functions on grpFrom
	for(size_t from = 0; from < grpFromSmall.size(); ++from)
	{
	//	get unique id of function
		const size_t uniqueID = grpFromSmall[from];

	//	find unique id of function in grpTo
		const size_t locIndex = grpToLarge.local_index(uniqueID);

	//	set mapping
		map.push_back(locIndex);
	}
}

void CreateFunctionIndexMapping(FunctionIndexMapping& map,
                                const FunctionGroup& grpFrom,
                                ConstSmartPtr<FunctionPattern> fctPattern)
{
	FunctionGroup commonFctGroup(fctPattern);
	commonFctGroup.add_all();
	CreateFunctionIndexMapping(map, grpFrom, commonFctGroup);
}


/**
 * This function create the union of function groups. Container is clear at beginning.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vFctGrp		Vector of function group (may contain nullptr)
 * \param[in]		sortFct		flag if group should be sorted after adding
 */
void CreateUnionOfFunctionGroups(FunctionGroup& fctGrp,
                                 const vector<const FunctionGroup*>& vFctGrp,
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
		if(vFctGrp[grp] == nullptr) continue;

		ConstSmartPtr<FunctionPattern> pFctPat = vFctGrp[grp]->function_pattern();
		if(pFctPat.invalid())
			UG_THROW("CreateUnionOfFunctionGroups: Function group "
					<<grp<<" has nullptr as underlying FunctionPattern.");

		fctGrp.set_function_pattern(pFctPat);
		break;
	}

//	if no function group given
	if(grp == vFctGrp.size()) return;

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vFctGrp.size(); ++i)
	{
	//	add subset group of elem disc
		if(vFctGrp[i] != nullptr)
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

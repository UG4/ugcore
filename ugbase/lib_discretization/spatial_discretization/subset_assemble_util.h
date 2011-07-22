/*
 * subset_assemble_util.h
 *
 *  Created on: 08.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__
#define __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_discretization/common/function_group.h"
#include "lib_discretization/common/subset_group.h"
#include "lib_discretization/common/groups_util.h"


namespace ug {

/**
 * This function create the union of subsets on which the Element Discretizations
 * work.
 *
 * \param[out]		ssGrp		Union of Subsets
 * \param[in]		vElemDisc	Vector of element discs, defined for subsets
 * \param[in]		clearGroup	flag if group should be cleared
 */
inline
bool CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp,
                        SubsetGroup& unionSSGrp,
                        std::vector<IElemDisc* > vElemDisc,
                        const ISubsetHandler& sh)
{
//	resize subset group vector
	vSSGrp.resize(vElemDisc.size());

//	if empty, nothing to do
	if(vSSGrp.empty()) {unionSSGrp.clear(); return true;}

//	create subset group for each elem disc
	for(size_t i = 0; i < vSSGrp.size(); ++i)
	{
	//	create subset group for elem disc i
		if(!ConvertStringToSubsetGroup(vSSGrp[i], sh,
		                           vElemDisc[i]->symb_subsets()))
		{
			UG_LOG("ERROR in 'CreateUnionOfSubsets': Cannot find symbolic "
					" subset name for IElemDisc "<<i<<".\n");
			return false;
		}
	}

//	set underlying subsetHandler
	unionSSGrp.set_subset_handler(sh);

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vSSGrp.size(); ++i)
	{
	//	add subset group of elem disc
		if(!unionSSGrp.add(vSSGrp[i]))
		{
			UG_LOG("ERROR in 'CreateUnionOfSubsets': Cannot add subsets of the "
				   "Elem Disc "<< i << " to union of Subsets.\n");
			return false;
		}
	}

//	we're done
	return true;
}

/**
 * This function selects from a given set of element discretizations those
 * who work on a given subset
 *
 * \param[out]		vSubsetElemDisc		Elem Disc working on subset
 * \param[in]		vElemDisc			Elem Disc to select from
 * \param[in]		si					Subset index
 * \param[in]		clearVec			flag if vector should be cleared
 */
inline
bool GetElemDiscOnSubset(std::vector<IElemDisc*>& vSubsetElemDisc,
                         const std::vector<IElemDisc*>& vElemDisc,
                         const std::vector<SubsetGroup>& vSSGrp,
                         int si, bool clearVec = true)
{
//	clear Vector
	if(clearVec) vSubsetElemDisc.clear();

//	loop elem discs
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
	//	if subset is used, add elem disc to Subset Elem Discs
		if(vSSGrp[i].contains(si))
			vSubsetElemDisc.push_back(vElemDisc[i]);
	}

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__ */

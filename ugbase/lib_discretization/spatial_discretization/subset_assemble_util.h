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


namespace ug {

/**
 * This function create the union of subsets on which the Element Discretizations
 * work.
 *
 * \param[out]		ssGrp		Union of Subsets
 * \param[in]		vElemDisc	Vector of element discs, defined for subsets
 * \param[in]		clearGroup	flag if group should be cleared
 */
template <typename TAlgebra>
bool CreateUnionOfSubsets(SubsetGroup& ssGrp,
                          const std::vector<IElemDisc<TAlgebra>*>& vElemDisc)
{
//	if empty, nothing to do
	if(vElemDisc.empty()) {ssGrp.clear(); return true;}

//	set underlying subsetHandler
	ssGrp.set_subset_handler(*vElemDisc[0]->get_subset_group().get_subset_handler());

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
	//	add subset group of elem disc
		if(!ssGrp.add(vElemDisc[i]->get_subset_group()))
		{
			UG_LOG("ERROR in 'CreateUnionOfSubsets': Cannot add subsets of the "
				   "Elem Disc "<< i << ".\n");
			return false;
		}
	}

//	we're done
	return true;
}

/**
 * This function create the union of functions on which the Element Discretizations
 * work.
 *
 * \param[out]		fctGrp		Union of Functions
 * \param[in]		vElemDisc	Vector of element discs
 * \param[in]		clearGroup	flag if group should be cleared
 * \param[in]		clearGroup	flag if group should be sorted after adding
 */
template <typename TAlgebra>
bool CreateUnionOfFunctions(FunctionGroup& fctGrp,
                            const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
                            bool sortFct = false)
{
//	if empty, nothing to do
	if(vElemDisc.empty()) {fctGrp.clear(); return true;}

//	set underlying subsetHandler
	fctGrp.set_function_pattern(*vElemDisc[0]->get_function_group().get_function_pattern());

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
	//	add subset group of elem disc
		if(!fctGrp.add(vElemDisc[i]->get_function_group()))
		{
			UG_LOG("ERROR in 'CreateUnionOfSubsets': Cannot add functions of the "
				   "Elem Disc "<< i << ".\n");
			return false;
		}
	}

//	sort iff required
	if(sortFct) fctGrp.sort();

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
template <typename TAlgebra>
bool GetElemDiscOnSubset(std::vector<IElemDisc<TAlgebra>*>& vSubsetElemDisc,
                         const std::vector<IElemDisc<TAlgebra>*>& vElemDisc,
                         int si, bool clearVec = true)
{
//	clear Vector
	if(clearVec) vSubsetElemDisc.clear();

//	loop elem discs
	for(size_t i = 0; i < vElemDisc.size(); ++i)
	{
	//	get subset group of elem disc
		const SubsetGroup& ssGrp = vElemDisc[i]->get_subset_group();

	//	if subset is used, add elem disc to Subset Elem Discs
		if(ssGrp.contains(si))
			vSubsetElemDisc.push_back(vElemDisc[i]);
	}

//	we're done
	return true;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__SPATIAL_DISCRETIZATION__SUBSET_ASSEMBLE_UTIL__ */

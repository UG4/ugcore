/*
 * subset_assemble_util.h
 *
 *  Created on: 08.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SUBSET_ASSEMBLE_UTIL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SUBSET_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug {

/**
 * This function create the union of subsets on which the Element Discretizations
 * work.
 *
 * \param[out]		ssGrp		Union of Subsets
 * \param[in]		vElemDisc	Vector of element discs, defined for subsets
 * \param[in]		clearGroup	flag if group should be cleared
 */
template <typename TDomain>
void CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp,
                        SubsetGroup& unionSSGrp,
                        std::vector<IElemDisc<TDomain>* > vElemDisc,
                        ConstSmartPtr<ISubsetHandler> pSH);

/**
 * This function selects from a given set of element discretizations those
 * who work on a given subset
 *
 * \param[out]		vSubsetElemDisc		Elem Disc working on subset
 * \param[in]		vElemDisc			Elem Disc to select from
 * \param[in]		si					Subset index
 * \param[in]		clearVec			flag if vector should be cleared
 */
template <typename TDomain>
void GetElemDiscOnSubset(std::vector<IElemDisc<TDomain>*>& vSubsetElemDisc,
                         const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                         const std::vector<SubsetGroup>& vSSGrp,
                         int si, bool clearVec = true);

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__SUBSET_ASSEMBLE_UTIL__ */

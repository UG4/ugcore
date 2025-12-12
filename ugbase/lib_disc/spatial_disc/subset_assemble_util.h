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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__SUBSET_ASSEMBLE_UTIL__
#define __H__UG__LIB_DISC__SPATIAL_DISC__SUBSET_ASSEMBLE_UTIL__

// extern includes
#include <iostream>
#include <vector>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/tools/subset_group.h"
// #include "lib_disc/common/function_group.h"
#include "lib_disc/common/groups_util.h"
#include "lib_disc/spatial_disc/elem_disc/elem_disc_interface.h"

namespace ug {

/**
 * This function create the union of subsets on which the Element Discretizations
 * work.
 *
 * \param[out]		vSSGrp		Union of Subsets
 * \param[in]		vElemDisc	Vector of element discs, defined for subsets
 * \param[in]		clearGroup	flag if group should be cleared
 */

template <typename TElemIn>
void CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp,
                        SubsetGroup& unionSSGrp,
                        std::vector<TElemIn*> vElemDisc,
                        ConstSmartPtr<ISubsetHandler> pSH)
{
	PROFILE_FUNC_GROUP("discretization");
//	resize subset group vector
	vSSGrp.resize(vElemDisc.size());

//	if empty, nothing to do
	if(vSSGrp.empty()) {unionSSGrp.clear(); return;}

//	create subset group for each elem disc
	for(size_t i = 0; i < vSSGrp.size(); ++i)
	{
		vSSGrp[i].set_subset_handler(pSH);
		vSSGrp[i].add(vElemDisc[i]->symb_subsets());
	}

//	set underlying subsetHandler
	unionSSGrp.set_subset_handler(pSH);

//	add all Subset groups of the element discs
	for(size_t i = 0; i < vSSGrp.size(); ++i)
	{
		//	add subset group of elem disc
		try{
			unionSSGrp.add(vSSGrp[i]);
		}UG_CATCH_THROW("Cannot add subsets of the Elem Disc "<<i<<" to union of Subsets.");
	}
}




/**
 * This function selects from a given set of element discretizations those
 * who work on a given subset.
 *
 *  At the same time, pointers are converted
 *  e.g. from TElemIn=IElemDisc<TDomain> to typename TElemOut=IElemError<TDomain>
 *
 * \param[out]		vSubsetElemDiscOut		Elem Disc items  working on subset (pointers type of type TElemOut)
 * \param[in]		vElemDiscIn			Elem Disc items to select from (pointers type of type TElemIn)
 * \param[in]		si					Subset index
 * \param[in]		clearVec			flag if vector should be cleared
 */
template <typename TElemOut, typename TElemIn>
void GetElemDiscItemOnSubset(std::vector<TElemOut*>& vSubsetElemDiscOut,
                         const std::vector<TElemIn*>& vElemDiscIn,
                         const std::vector<SubsetGroup>& vSSGrp,
                         int si, bool clearVec=true)
{
//	clear Vector
	if(clearVec) vSubsetElemDiscOut.clear();

//	loop elem discs
	for(size_t i = 0; i < vElemDiscIn.size(); ++i)
	{
	//	if subset is used, add elem disc to Subset Elem Discs
		if(vSSGrp[i].contains(si))
			vSubsetElemDiscOut.push_back(static_cast<TElemOut*>(vElemDiscIn[i]));
	}
}

/**
 * Legacy version of 'GetElemDiscItemOnSubset' for TElemOut=TElemIn=IElemDisc<TDomain>
 */
template <typename TDomain>
void GetElemDiscOnSubset(std::vector<IElemDisc<TDomain>*>& vSubsetElemDisc,
                         const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                         const std::vector<SubsetGroup>& vSSGrp,
                         int si, bool clearVec = true);


} // end namespace ug

#endif

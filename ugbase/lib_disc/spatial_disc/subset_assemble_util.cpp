/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
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

#include "subset_assemble_util.h"

namespace ug{

template <typename TDomain>
void CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp,
                        SubsetGroup& unionSSGrp,
                        std::vector<IElemDisc<TDomain>* > vElemDisc,
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

template <typename TDomain>
void GetElemDiscOnSubset(std::vector<IElemDisc<TDomain>*>& vSubsetElemDisc,
                         const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                         const std::vector<SubsetGroup>& vSSGrp,
                         int si, bool clearVec)
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
}

////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template void CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp, SubsetGroup& unionSSGrp, std::vector<IElemDisc<Domain1d>* > vElemDisc, ConstSmartPtr<ISubsetHandler> pSH);
template void GetElemDiscOnSubset(std::vector<IElemDisc<Domain1d>*>& vSubsetElemDisc, const std::vector<IElemDisc<Domain1d>*>& vElemDisc, const std::vector<SubsetGroup>& vSSGrp, int si, bool clearVec);
#endif
#ifdef UG_DIM_2
template void CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp, SubsetGroup& unionSSGrp, std::vector<IElemDisc<Domain2d>* > vElemDisc, ConstSmartPtr<ISubsetHandler> pSH);
template void GetElemDiscOnSubset(std::vector<IElemDisc<Domain2d>*>& vSubsetElemDisc, const std::vector<IElemDisc<Domain2d>*>& vElemDisc, const std::vector<SubsetGroup>& vSSGrp, int si, bool clearVec);
#endif
#ifdef UG_DIM_3
template void CreateSubsetGroups(std::vector<SubsetGroup>& vSSGrp, SubsetGroup& unionSSGrp, std::vector<IElemDisc<Domain3d>* > vElemDisc, ConstSmartPtr<ISubsetHandler> pSH);
template void GetElemDiscOnSubset(std::vector<IElemDisc<Domain3d>*>& vSubsetElemDisc, const std::vector<IElemDisc<Domain3d>*>& vElemDisc, const std::vector<SubsetGroup>& vSSGrp, int si, bool clearVec);
#endif

} // end namespace ug

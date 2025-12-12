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

namespace ug {



template <typename TDomain>
void GetElemDiscOnSubset(std::vector<IElemDisc<TDomain>*>& vSubsetElemDisc,
                         const std::vector<IElemDisc<TDomain>*>& vElemDisc,
                         const std::vector<SubsetGroup>& vSSGrp,
                         int si, bool clearVec)
{
	GetElemDiscItemOnSubset<IElemDisc<TDomain>, IElemDisc<TDomain> >
	(vSubsetElemDisc, vElemDisc, vSSGrp, si, clearVec);
}



////////////////////////////////////////////////////////////////////////////////
//	explicit template instantiations
////////////////////////////////////////////////////////////////////////////////

#ifdef UG_DIM_1
template void CreateSubsetGroups<IElemDisc<Domain1d> >(std::vector<SubsetGroup>& vSSGrp, SubsetGroup& unionSSGrp, std::vector<IElemDisc<Domain1d>* > vElemDisc, ConstSmartPtr<ISubsetHandler> pSH);
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

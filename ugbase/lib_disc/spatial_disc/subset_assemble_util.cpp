
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

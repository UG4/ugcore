/*
 * dof_info.cpp
 *
 *  Created on: 15.02.2013
 *      Author: andreasvogel
 */

#include "dof_distribution_info.h"
#include "lib_disc/common/function_group.h"
#include "lib_disc/common/subset_group.h"
#include "lib_disc/reference_element/reference_element_util.h"
#include "lib_disc/local_finite_element/local_finite_element_provider.h"

#include "common/log.h"
#include <iostream>
using namespace std;

namespace ug{

DoFDistributionInfo::DoFDistributionInfo(ConstSmartPtr<ISubsetHandler> spSH)
	: FunctionPattern(spSH)
{}

void DoFDistributionInfo::init()
{
	PROFILE_FUNC();
	FunctionPattern::lock();
	create_offsets();
}

void DoFDistributionInfo::create_offsets(ReferenceObjectID roid)
{
	PROFILE_FUNC();
// 	loop subsets
	for(int si = 0; si < num_subsets(); ++si)
	{
	//	reset
		m_vvNumDoFsOnROID[roid][si] = 0;

	//	loop functions
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
		//	if function is not defined, we leave the offset as invalid.
			if(!is_def_in_subset(fct, si)) {
				m_vvvOffsets[roid][si][fct] = NOT_DEF_ON_SUBSET;
				continue;
			}

		//	get trial space
			const CommonLocalDoFSet& clds = LocalFiniteElementProvider::get_dofs(lfeid(fct));

		//	get number of DoFs on the reference element need for the space
			const int numDoF = clds.num_dof(roid);

		//	check that numDoFs specified by this roid
			if(numDoF == NOT_SPECIFIED || m_vvNumDoFsOnROID[roid][si]){
				m_vvNumDoFsOnROID[roid][si] = NOT_SPECIFIED;
				m_vvvOffsets[roid][si][fct] = 0;
				continue;
			}

		//	overwrite max dim with dofs (if subset has that dimension)
			if(numDoF > 0)
				if(dim_subset(si) > ReferenceElementDimension(roid))
					m_vMaxDimToOrderDoFs[fct] = ReferenceElementDimension(roid);

		//	set offset for each function defined in the subset
			m_vvvOffsets[roid][si][fct] = m_vvNumDoFsOnROID[roid][si];

		//	increase number of DoFs per type on this subset
			m_vvNumDoFsOnROID[roid][si] += numDoF;
		}
	}
}

void DoFDistributionInfo::create_offsets()
{
	PROFILE_FUNC();
//	function infos
	m_vMaxDimToOrderDoFs.resize(num_fct(), 0);

//	cache number of DoFs in a sub-geometric object
	for(size_t fct = 0; fct < num_fct(); ++fct){
		const CommonLocalDoFSet& lds = LocalFiniteElementProvider::get_dofs(lfeid(fct));

		for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid){
			m_vNumDoFOnSubelem[roid].resize(num_fct());
			m_vNumDoFOnSubelem[roid][fct] = lds.num_dof((ReferenceObjectID)roid);
		}
	}

//	loop all reference element to resize the arrays
	for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid)
	{
	//	clear offsets
		m_vvvOffsets[roid].clear();
		m_vvNumDoFsOnROID[roid].clear();

	//	resize for all subsets
		m_vvvOffsets[roid].resize(num_subsets());
		m_vvNumDoFsOnROID[roid].resize(num_subsets(), 0);

	//	resize for each function
		for(int si = 0; si < num_subsets(); ++si)
			m_vvvOffsets[roid][si].resize(num_fct());
	}

//	loop all reference element, but not vertices (no disc there)
	for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid)
		create_offsets((ReferenceObjectID) roid);

//	reset dimension maximum
	for(size_t d=0; d < NUM_GEOMETRIC_BASE_OBJECTS; ++d)
	{
		m_vMaxDoFsInDim[d] = 0;
		m_vvMaxDoFsInDimPerSubset[d].resize(num_subsets(), 0);
	}

//	get max number of dofs per roid
	for(int roid=ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid)
	{
	//	lets find out the maximum number of DoFs for objects in dimension
		const int d = ReferenceElementDimension((ReferenceObjectID)roid);

	// 	lets find out the maximum number for a type on all subsets
		m_vMaxDoFsOnROID[roid] = 0;
		for(int si = 0; si < num_subsets(); ++si)
		{
			if(m_vvNumDoFsOnROID[roid][si] != NOT_SPECIFIED)
				m_vMaxDoFsOnROID[roid] = std::max(m_vMaxDoFsOnROID[roid],
												  (size_t)m_vvNumDoFsOnROID[roid][si]);

			if(m_vvNumDoFsOnROID[roid][si] != NOT_SPECIFIED)
				m_vvMaxDoFsInDimPerSubset[d][si] = std::max(m_vvMaxDoFsInDimPerSubset[d][si],
				                                            (size_t)m_vvNumDoFsOnROID[roid][si]);
		}

	//	compute maximum per dim objects
		m_vMaxDoFsInDim[d] = std::max(m_vMaxDoFsInDim[d],
		                              m_vMaxDoFsOnROID[roid]);
	}
}

void DoFDistributionInfo::print_local_dof_statistic(int verboseLev) const
{
//	Subset informations
	UG_LOG(num_subsets() << " Subset(s) used (Subset Name, dim): ");
	for(int si = 0; si < num_subsets(); ++si)
	{
		if(si > 0) UG_LOG(", ");
		UG_LOG("(" << subset_name(si) << ", " << dim_subset(si) << ")");
	}
	UG_LOG("\n");

//	Function informations
	UG_LOG(num_fct() << " Function(s) defined (Symbolic Name, dim): ");
	for(size_t fct = 0; fct < num_fct(); ++fct)
	{
		if(fct > 0) UG_LOG(", ");
		UG_LOG("(" << name(fct) << ", " << dim(fct) << ")");
	}
	UG_LOG("\n");

//	print subsets of functions
	if(verboseLev >= 2)
	{
		UG_LOG("Function definition on subsets: \n");
		for(size_t fct = 0; fct < num_fct(); ++fct)
		{
			UG_LOG("   "<<name(fct) << ": ");
			if(is_def_everywhere(fct)) UG_LOG("[everywhere] ");
			bool bFirst = true;
			for(int si = 0; si < num_subsets(); ++si)
			{
				if(bFirst) bFirst = false; else UG_LOG(", ");
				if(!is_def_in_subset(fct, si)) continue;
				UG_LOG(subset_name(si));
			}
			UG_LOG("\n");
		}
		UG_LOG("\n");
	}

//	write info for subset/fct -> localFEId info
	UG_LOG("\n\t\t\t Subsets\n");
	UG_LOG(" "<<setw(14)<<"Function"<<" |");
	for(int si = 0; si < num_subsets(); ++si)
		UG_LOG(setw(11)<<si<<setw(8)<<" "<<"|")
	UG_LOG("\n");
	for(size_t fct = 0; fct < num_fct(); ++fct)
	{
		UG_LOG(" "<<setw(14)<<name(fct)<<" |");
		for(int si = 0; si < num_subsets(); ++si)
		{
			if(!is_def_in_subset(fct,si))
				 {UG_LOG(setw(8)<<"---"<<setw(8)<<" "<<"|");}
			else {UG_LOG(setw(16)<<lfeid(fct)<<" |");}
		}
		UG_LOG("\n");
	}

//	write info about DoFs on ROID
	UG_LOG("\n");
	UG_LOG("                  | "<<"        "<<"  |  Subsets \n");
	UG_LOG(" ReferenceElement |");
	UG_LOG("   "<<setw(4)<<"max"<<"    |");
	for(int si = 0; si < num_subsets(); ++si)
		UG_LOG("   "<<setw(4)<<si<<"    |");
	UG_LOG("\n")
	UG_LOG("-------------------");
	for(int si = 0; si <= num_subsets(); ++si)
		UG_LOG("-------------");
	UG_LOG("\n")

	for(int i=ROID_VERTEX; i < NUM_REFERENCE_OBJECTS; ++i)
	{
		ReferenceObjectID roid = (ReferenceObjectID) i;

		UG_LOG(" " << setw(16) << roid << " |");
		UG_LOG("   "<<setw(4) << m_vMaxDoFsOnROID[roid] << "    |");
		for(int si = 0; si < num_subsets(); ++si)
			UG_LOG("   "<<setw(4) << m_vvNumDoFsOnROID[roid][si] << "    |");

		UG_LOG("\n");
	}
	for(int d = 0; d < NUM_GEOMETRIC_BASE_OBJECTS; ++d)
	{
		UG_LOG(setw(14) << " all " <<setw(2)<< d << "d |");
		UG_LOG("   "<<setw(4) << m_vMaxDoFsInDim[d] << "    |");
		UG_LOG("\n");
	}
	UG_LOG("\n");
}


FunctionGroup DoFDistributionInfoProvider::fct_grp_by_name(const char* names) const
{
	return FunctionGroup(m_spDDI, TokenizeString(names));
}

SubsetGroup DoFDistributionInfoProvider::subset_grp_by_name(const char* names) const
{
	return m_spDDI->subset_grp_by_name(names);
}


} // end namespace ug

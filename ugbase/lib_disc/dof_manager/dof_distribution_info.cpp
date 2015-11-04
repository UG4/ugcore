
#include "dof_distribution_info.h"
#include "lib_grid/tools/subset_group.h"
#include "lib_disc/common/function_group.h"
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


void DoFDistributionInfo::create_offsets()
{
	PROFILE_FUNC();

	/////////////////////////////////////////////
	// resize arrays and reset values to zero
	/////////////////////////////////////////////

//	resize for subsets
	for(int gbo = VERTEX; gbo < NUM_GEOMETRIC_BASE_OBJECTS; ++gbo){
		m_vMaxDoFsInDim[gbo] = 0;
		m_vvMaxDoFsInDimPerSubset[gbo].clear();
		m_vvMaxDoFsInDimPerSubset[gbo].resize(num_subsets(), 0);
	}

	for(int roid = ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid){
		m_vMaxDoFsOnROID[roid] = 0;
		m_vvNumDoFsOnROIDPerSubset[roid].clear();
		m_vvNumDoFsOnROIDPerSubset[roid].resize(num_subsets(), 0);
	}

//	resize function infos
	m_vFctInfo.clear();
	m_vFctInfo.resize(num_fct());

	for(size_t fct = 0; fct < num_fct(); ++fct){

		for(int gbo = VERTEX; gbo < NUM_GEOMETRIC_BASE_OBJECTS; ++gbo){
			m_vFctInfo[fct].vMaxDoFsInDim[gbo] = 0;
			m_vFctInfo[fct].vvMaxDoFsInDimPerSubset[gbo].clear();
			m_vFctInfo[fct].vvMaxDoFsInDimPerSubset[gbo].resize(num_subsets(), 0);
		}

		for(int roid = ROID_VERTEX; roid < NUM_REFERENCE_OBJECTS; ++roid){
			m_vFctInfo[fct].vMaxDoFsOnROID[roid] = 0;
			m_vFctInfo[fct].vvNumDoFsOnROIDPerSubset[roid].clear();
			m_vFctInfo[fct].vvNumDoFsOnROIDPerSubset[roid].resize(num_subsets(), 0);

			m_vFctInfo[fct].vvOffsets[roid].resize(num_subsets(), NOT_DEF_ON_SUBSET);
		}
	}

	/////////////////////////////////////////////
	// compute values
	/////////////////////////////////////////////

//	loop reference element by reference element
	for(int r = ROID_VERTEX; r < NUM_REFERENCE_OBJECTS; ++r)
	{
	//	get reference element and dimension
		const ReferenceObjectID roid = (ReferenceObjectID) r;
		const int d = ReferenceElementDimension(roid);

	//	loop subsets and functions
		for(int si = 0; si < num_subsets(); ++si)
		{
			for(size_t fct = 0; fct < num_fct(); ++fct)
			{
			//	get dof info for space
				const CommonLocalDoFSet& lds =
							LocalFiniteElementProvider::get_dofs(lfeid(fct));

			//	function info
				FctInfo& fi = m_vFctInfo[fct];

			//	get numDoFs on (fct,roid,si)
				size_t numFctDoFOnROIDPerSubset = 0;

			//	if number specified by the space, set it
				if(lds.num_dof(roid) != CommonLocalDoFSet::NOT_SPECIFIED)
					numFctDoFOnROIDPerSubset = (size_t)lds.num_dof(roid);

			//	if function not defined on subset, number of DoFs is zero
				if(!is_def_in_subset(fct, si))
					numFctDoFOnROIDPerSubset = 0;


				///// Infos on fct //////


				// set num DoFs (always >= 0)
				fi.vvNumDoFsOnROIDPerSubset[roid][si] = numFctDoFOnROIDPerSubset;

				// get fct maximum on roid on whole grid (always >= 0)
				fi.vMaxDoFsOnROID[roid] =
					std::max(fi.vMaxDoFsOnROID[roid], numFctDoFOnROIDPerSubset);

				// get fct maximum in dimension per subset (always >= 0)
				fi.vvMaxDoFsInDimPerSubset[d][si] =
					std::max(fi.vvMaxDoFsInDimPerSubset[d][si], numFctDoFOnROIDPerSubset);

				// get fct maximum in dimension on whole grid (always >= 0)
				fi.vMaxDoFsInDim[d] =
					std::max(fi.vMaxDoFsInDim[d], numFctDoFOnROIDPerSubset);


				///// Summarize over functions //////

				//	set offset for each function defined in the subset
				fi.vvOffsets[roid][si] = m_vvNumDoFsOnROIDPerSubset[roid][si];

				// count number of dofs on (roid, si)
				m_vvNumDoFsOnROIDPerSubset[roid][si] += numFctDoFOnROIDPerSubset;

			} // end fct
		} // end subset


		// loop subsets and sum up global maxima
		for(int si = 0; si < num_subsets(); ++si)
		{
			const size_t numDoFOnROIDPerSubset = m_vvNumDoFsOnROIDPerSubset[roid][si];

			// get maximum in dimension on whole grid (always >= 0)
			m_vMaxDoFsOnROID[roid] =
				std::max(m_vMaxDoFsOnROID[roid], numDoFOnROIDPerSubset);

			// get maximum in dimension per subset (always >= 0)
			m_vvMaxDoFsInDimPerSubset[d][si] =
				std::max(m_vvMaxDoFsInDimPerSubset[d][si], numDoFOnROIDPerSubset);

			// get maximum on roid on whole grid (always >= 0)
			m_vMaxDoFsInDim[d] =
				std::max(m_vMaxDoFsInDim[d], numDoFOnROIDPerSubset);
		} // end subset

	} // end roid
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
			UG_LOG("   "<<setw(4) << m_vvNumDoFsOnROIDPerSubset[roid][si] << "    |");

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

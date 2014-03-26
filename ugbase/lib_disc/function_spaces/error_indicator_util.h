
#ifndef __H__UG_DISC__ERROR_INDICATOR_UTIL__
#define __H__UG_DISC__ERROR_INDICATOR_UTIL__

#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug{

/// marks elements according to an attached error value field
/**
 * This function marks elements for refinement. The passed error attachment
 * is used as a weight for the amount of the error an each element. All elements
 * that have an indicated error with s* max <= err <= max are marked for refinement.
 * Here, max is the maximum error measured, s is a scaling quantity chosen by
 * the user. In addition, all elements with an error smaller than TOL
 * (user defined) are not refined.
 *
 * \param[in, out]	refiner		Refiner, elements marked on exit
 * \param[in]		dd			dof distribution
 * \param[in]		TOL			Minimum error, such that an element is marked
 * \param[in]		scale		scaling factor indicating lower bound for marking
 * \param[in]		aaError		Error value attachment to elements
 */
template <typename TElem>
void MarkElements(MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
                  IRefiner& refiner,
                  ConstSmartPtr<DoFDistribution> dd,
                  number TOL,
                  number refineFrac, number coarseFrac,
				  int maxLevel)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

//	reset maximum of error
	number max = 0.0, min = std::numeric_limits<number>::max();
	number totalErr = 0.0;
	int numElem = 0;

//	get element iterator
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop all elements to find the maximum of the error
	for( ;iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	search for maximum and minimum
		if(aaError[elem] > max)	max = aaError[elem];
		if(aaError[elem] < min)	min = aaError[elem];

	//	sum up total error
		totalErr += aaError[elem];
		numElem += 1;
	}

#ifdef UG_PARALLEL
	number maxLocal = max, minLocal = min, totalErrLocal = totalErr;
	int numElemLocal = numElem;
	if(pcl::NumProcs() > 1){
		pcl::ProcessCommunicator com;
		max = com.allreduce(maxLocal, PCL_RO_MAX);
		min = com.allreduce(minLocal, PCL_RO_MIN);
		totalErr = com.allreduce(totalErrLocal, PCL_RO_SUM);
		numElem = com.allreduce(numElemLocal, PCL_RO_SUM);
	}
#endif
	UG_LOG("  +++++  Gradient Error Indicator on "<<numElem<<" Elements +++++\n");
#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG("  +++ Element Errors on Proc " << pcl::ProcRank() <<
			   ": maximum=" << max << ", minimum="<<min<<", sum="<<totalErr<<".\n");
	}
#endif

	UG_LOG("  +++ Element Errors: maximum=" << max << ", minimum="<<min<<
	       ", sum="<<totalErr<<".\n");

//	check if total error is smaller than tolerance. If that is the case we're done
	if(totalErr < TOL)
	{
		UG_LOG("  +++ Total error "<<totalErr<<" smaller than TOL ("<<TOL<<"). done.");
		return;
	}

//	Compute minimum
	number minErrToRefine = max * refineFrac;
	UG_LOG("  +++ Refining elements if error greater " << refineFrac << "*" <<max<<
			" = "<< minErrToRefine << ".\n");
	number maxErrToCoarse = min * (1+coarseFrac);
	if(maxErrToCoarse < TOL/numElem) maxErrToCoarse = TOL/numElem;
	UG_LOG("  +++ Coarsening elements if error smaller "<< maxErrToCoarse << ".\n");

//	reset counter
	int numMarkedRefine = 0, numMarkedCoarse = 0;

	iter = dd->template begin<TElem>();
	iterEnd = dd->template end<TElem>();

//	loop elements for marking
	for(; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	marks for refinement
		if(aaError[elem] >= minErrToRefine)
			if(dd->multi_grid()->get_level(elem) <= maxLevel)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}

	//	marks for coarsening
		if(aaError[elem] <= maxErrToCoarse)
		{
			refiner.mark(elem, RM_COARSEN);
			numMarkedCoarse++;
		}
	}

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG("  +++ Marked for refinement on Proc "<<pcl::ProcRank()<<": " << numMarkedRefine << " Elements.\n");
		UG_LOG("  +++ Marked for coarsening on Proc "<<pcl::ProcRank()<<": " << numMarkedCoarse << " Elements.\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " Elements.\n");
	UG_LOG("  +++ Marked for coarsening: " << numMarkedCoarse << " Elements.\n");
}

/// marks elements according to an attached error value field
/**
 * This function marks elements for refinement. The passed error attachment
 * is used as a weight for the amount of the error an each element. All elements
 * that have an indicated error > refineTol are marked for refinement and
 * elements with an error < coarsenTol are marked for coarsening
 *
 * \param[in, out]	refiner		Refiner, elements marked on exit
 * \param[in]		dd			dof distribution
 * \param[in]		refTol		all elements with error > refTol are marked for refinement.
 * 								If refTol is negative, no element will be marked for refinement.
 * \param[in]		coarsenTol	all elements with error < coarsenTol are marked for coarsening.
 * 								If coarsenTol is negative, no element will be marked for coarsening.
 * \param[in]		aaError		Error value attachment to elements
 */
template <typename TElem>
void MarkElementsAbsolute(MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
						  IRefiner& refiner,
						  ConstSmartPtr<DoFDistribution> dd,
						  number refTol,
						  number coarsenTol,
						  int minLevel,
						  int maxLevel)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	int numMarkedRefine = 0, numMarkedCoarse = 0;
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop elements for marking
	for(; iter != iterEnd; ++iter){
		TElem* elem = *iter;

	//	marks for refinement
		if((refTol >= 0)
			&& (aaError[elem] > refTol)
			&& (dd->multi_grid()->get_level(elem) < maxLevel))
		{
			refiner.mark(elem, RM_REFINE);
			numMarkedRefine++;
		}

	//	marks for coarsening
		if((coarsenTol >= 0)
			&& (aaError[elem] < coarsenTol)
			&& (dd->multi_grid()->get_level(elem) > minLevel))
		{
			refiner.mark(elem, RM_COARSEN);
			numMarkedCoarse++;
		}
	}

#ifdef UG_PARALLEL
	if(pcl::NumProcs() > 1){
		UG_LOG("  +++ Marked for refinement on Proc "<<pcl::ProcRank()<<": " << numMarkedRefine << " Elements.\n");
		UG_LOG("  +++ Marked for coarsening on Proc "<<pcl::ProcRank()<<": " << numMarkedCoarse << " Elements.\n");
		pcl::ProcessCommunicator com;
		int numMarkedRefineLocal = numMarkedRefine, numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " Elements.\n");
	UG_LOG("  +++ Marked for coarsening: " << numMarkedCoarse << " Elements.\n");
}

}//	end of namespace

#endif

/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Sebastian Reiter
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

#ifndef __H__UG_DISC__ERROR_INDICATOR_UTIL__
#define __H__UG_DISC__ERROR_INDICATOR_UTIL__

#include "lib_grid/multi_grid.h"
#include "lib_grid/refinement/refiner_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"

namespace ug{


/// helper function that computes min/max and total of error indicators
/**
 * @param[out] min		minimal eta_i
 * @param[out] max		maximal eta_i
 * @param[out] sum		sum of eta_i
 * @param[out] errSq	sum of eta_i^2
 * @param[out] numElem	number of elements
 * @return				average eta
 */
template<typename TElem>
number ComputeAvg
(	MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError2,
	ConstSmartPtr<DoFDistribution> dd,
	number& min, number& max, number& sum, number& errSq, size_t& numElem,
	number& minLocal, number& maxLocal, number& sumLocal, size_t& numElemLocal
)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

//	reset maximum of error
	max = 0.0, min = std::numeric_limits<number>::max();
	sum = 0.0;
	errSq = 0.0;
	numElem = 0;

//	get element iterator
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop all elements to find the maximum of the error
	for (; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

		const number elemEta2 = aaError2[elem];

	//	if no error value exists: ignore (might be newly added by refinement);
	//	newly added elements are supposed to have a negative error estimator
		if (aaError2[elem] < 0) UG_THROW("error value invalid!");//continue;


		const number elemEta = sqrt(aaError2[elem]);

	//	search for maximum and minimum
		if (elemEta > max) max = elemEta;
		if (elemEta < min) min = elemEta;

	//	sum up total error
		sum += elemEta;
		errSq += elemEta2;
		++numElem;
	}

	// set local variables
	maxLocal = max;
	minLocal = min;
	sumLocal = sum;
#ifdef UG_PARALLEL
	number errLocal = errSq;
#endif
	numElemLocal = numElem;

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		max = com.allreduce(maxLocal, PCL_RO_MAX);
		min = com.allreduce(minLocal, PCL_RO_MIN);
		sum = com.allreduce(sumLocal, PCL_RO_SUM);
		errSq = com.allreduce(errLocal, PCL_RO_SUM);
		numElem = com.allreduce(numElemLocal, PCL_RO_SUM);
	}
#endif
	UG_LOG("  +++++  Error indicator on " << numElem << " elements +++++\n");
	UG_LOG("  +++ Element errors: maxEta=" << max << ", minEta="
			<< min << ", sumEta=" << sum << ", sumEtaSq=" << errSq << ".\n");

	return (sum/numElem);
}


/// helper function that computes min/max and total of error indicators
/**
 * @param[out] min		minimal eta_i^2
 * @param[out] max		maximal eta_i^2
 * @param[out] totalErr	sum of eta_i^2
 * @param[out] numElem	number of elements
 *
 */
template<typename TElem>
void ComputeMinMax
(	MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
	ConstSmartPtr<DoFDistribution> dd,
	number& min, number& max, number& totalErr, size_t& numElem,
	number& minLocal, number& maxLocal, number& totalErrLocal, size_t& numElemLocal
)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

//	reset maximum of error
	max = 0.0, min = std::numeric_limits<number>::max();
	totalErr = 0.0;
	numElem = 0;

//	get element iterator
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop all elements to find the maximum of the error
	for (; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	if no error value exists: ignore (might be newly added by refinement);
	//	newly added elements are supposed to have a negative error estimator
		if (aaError[elem] < 0) continue;

	//	search for maximum and minimum
		if (aaError[elem] > max) max = aaError[elem];
		if (aaError[elem] < min) min = aaError[elem];

	//	sum up total error
		totalErr += aaError[elem];
		++numElem;
	}

	// set local variables
	maxLocal = max;
	minLocal = min;
	totalErrLocal = totalErr;
	numElemLocal = numElem;

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		max = com.allreduce(maxLocal, PCL_RO_MAX);
		min = com.allreduce(minLocal, PCL_RO_MIN);
		totalErr = com.allreduce(totalErrLocal, PCL_RO_SUM);
		numElem = com.allreduce(numElemLocal, PCL_RO_SUM);
	}
#endif
	UG_LOG("  +++++  Error indicator on " << numElem << " elements +++++\n");
	UG_LOG("  +++ Element errors: maxEtaSq=" << max << ", minEtaSq="
			<< min << ", sumEtaSq=" << totalErr << ".\n");
}


/// helper function that computes min/max and total of error indicators
template<typename TElem>
void ComputeMinMaxTotal
(	MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError2,
	ConstSmartPtr<DoFDistribution> dd,
	number& min, number& max, number& totalErr, size_t& numElem

)
{
	number minLocal, maxLocal, totalErrLocal;
	size_t numElemLocal;
	ComputeMinMax(aaError2, dd, min, max, totalErr, numElem, minLocal, maxLocal, totalErrLocal, numElemLocal);
}
/// marks elements according to an attached error value field
/**
 * This function marks elements for refinement and coarsening. The passed error attachment
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
 * \param[in]		aaError		Error value attachment to elements (\f$ \eta_i^2 \f$)
 */
template<typename TElem>
void MarkElements(MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
		IRefiner& refiner,
		ConstSmartPtr<DoFDistribution> dd,
		number TOL,
		number refineFrac, number coarseFrac,
		int maxLevel)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

// compute minimal/maximal/ total error and number of elements
	number min, max, totalErr;
	size_t numElem;
	ComputeMinMaxTotal(aaError, dd, min, max, totalErr, numElem);

//	check if total error is smaller than tolerance. If that is the case we're done
	if(totalErr < TOL)
	{
		UG_LOG("  +++ Total error "<<totalErr<<" smaller than TOL ("<<TOL<<"). done.");
		return;
	}

//	Compute minimum
	number minErrToRefine = max * refineFrac;
	UG_LOG("  +++ Refining elements if error greater " << refineFrac << "*" << max <<
			" = " << minErrToRefine << ".\n");
	number maxErrToCoarse = min * (1 + coarseFrac);
	if(maxErrToCoarse < TOL/numElem) maxErrToCoarse = TOL/numElem;
	UG_LOG("  +++ Coarsening elements if error smaller " << maxErrToCoarse << ".\n");

//	reset counter
	int numMarkedRefine = 0, numMarkedCoarse = 0;

	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

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
	if(pcl::NumProcs() > 1)
	{
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

/// marks elements according for refinement to an attached error value field
/**
 * This function marks elements for refinement. The passed error attachment
 * is used as an indicator for the the error on each element.
 * Elements are refined if the error sum taken over all elements is greater
 * than the tolerance value tol. In that case, elements with an indicated
 * error of err >= tol / #elems are marked for refinement if and only if their
 * multigrid level is below the tolerated maximum of maxLevel.
 *
 * \param[in]		aaError		error value attachment to elements (\f$ \eta_i^2 \f$)
 * \param[in, out]	refiner		refiner, elements marked on exit
 * \param[in]		dd			dof distribution
 * \param[in]		tol			tolerated global error (no refinement if error below)
 * \param[in]		maxLevel	maximal refinement level in multigrid
 */

template<typename TElem>
void MarkElementsForRefinement
(	MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
	IRefiner& refiner,
	ConstSmartPtr<DoFDistribution> dd,
	number tol,
	int maxLevel
)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

// compute minimal/maximal/ total error and number of elements
	number min, max, totalErr;
	size_t numElem;
	ComputeMinMaxTotal(aaError, dd, min, max, totalErr, numElem);

//	check if total error is smaller than tolerance; if that is the case we're done
	if (totalErr < tol)
	{
		UG_LOG("  +++ Total error "<<totalErr<<" smaller than TOL (" << tol << "). "
			   "No refinement necessary." << std::endl);
		return;
	}

//	compute minimum
	//number minErrToRefine = max * refineFrac;
	number minErrToRefine = tol / numElem;
	UG_LOG("  +++ Refining elements if error greater " << tol << "/" << numElem <<
		   " = " << minErrToRefine << ".\n");

//	reset counter
	size_t numMarkedRefine = 0;

	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop elements for marking
	for (; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	//	if no error value exists: ignore (might be newly added by refinement);
	//	newly added elements are supposed to have a negative error estimator
		if (aaError[elem] < 0) continue;

	//	marks for refinement
		if (aaError[elem] >= minErrToRefine)
			if (dd->multi_grid()->get_level(elem) < maxLevel)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		size_t numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ Marked for refinement: " << numMarkedRefine << " elements.\n");
}

/// marks elements for coarsening according to an attached error value field
/**
 * This function marks elements for coarsening. The passed error attachment
 * is used as an indicator for the the error on each element.
 * Elements one level below surface are marked if the average error of its
 * children is lower than tol / #elems / 4 / safety. The safety factor is
 * supposed to ensure that elements are not refined and coarsened back and
 * forth in a dynamic adaptive simulation.
 *
 * \param[in]		aaError		error value attachment to elements (\f$ \eta_i^2 \f$)
 * \param[in, out]	refiner		refiner, elements marked on exit
 * \param[in]		dd			dof distribution
 * \param[in]		tol			tolerated global error
 * \param[in]		safety		safety factor
 */
template<typename TElem>
void MarkElementsForCoarsening
(	MultiGrid::AttachmentAccessor<TElem,
	ug::Attachment<number> >& aaError,
		IRefiner& refiner,
		ConstSmartPtr<DoFDistribution> dd,
		number TOL,
		number safety,
		int minLevel = 0
)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

// compute minimal/maximal/ total error and number of elements
	number min, max, totalErr;
	size_t numElem;
	ComputeMinMaxTotal(aaError, dd, min, max, totalErr, numElem);

//	compute maximum
	//number maxErrToCoarse = min * (1+coarseFrac);
	//if (maxErrToCoarse < TOL/numElem) maxErrToCoarse = TOL/numElem;
	number maxErrToCoarse = TOL / numElem / 4.0 / safety; // = tol_per_elem / (2^2*safety_factor)
	UG_LOG("  +++ Coarsening elements if avg child error smaller than "<< maxErrToCoarse << ".\n");

//	reset counter
	size_t numMarkedCoarse = 0;

	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();

//	loop elements for marking
	for (; iter != iterEnd; ++iter)
	{
	//	get element
		TElem* elem = *iter;

	// ignore if already marked for coarsening
		if (refiner.get_mark(elem) & RM_COARSEN) continue;

	// ignore if level too low
		if (dd->multi_grid()->get_level(elem) <= minLevel)
			 continue;

	// get parent
		TElem* parent = dynamic_cast<TElem*>(dd->multi_grid()->get_parent(elem));
		if (!parent) continue; // either dynamic casting failed or parent does not exist

	// sum up error values over all children of parent
		number sum = 0.0;
		size_t cnt = 0;
		size_t nCh = dd->multi_grid()->num_children<TElem>(parent);
		for (size_t ch = 0; ch < nCh; ch++)
		{
			TElem* child = dd->multi_grid()->get_child<TElem>(parent, ch);

			//	if no error value exists: ignore (might be newly added by refinement);
			//	newly added elements are supposed to have a negative error estimator
			if (aaError[child] < 0) continue;

			sum += aaError[child];
			cnt++;
		}

	//	marks for coarsening
		if (sum/cnt <= maxErrToCoarse)
		{
			for (size_t ch = 0; ch < nCh; ch++)
			{
				TElem* child = dd->multi_grid()->get_child<TElem>(parent, ch);
				if (refiner.get_mark(child) & RM_COARSEN) continue;
				refiner.mark(child, RM_COARSEN);
				numMarkedCoarse++;
			}
		}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		size_t numMarkedCoarseLocal = numMarkedCoarse;
		numMarkedCoarse = com.allreduce(numMarkedCoarseLocal, PCL_RO_SUM);
	}
#endif

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
template<typename TElem>
void MarkElementsAbsolute(MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> >& aaError,
						  IRefiner& refiner,
						  ConstSmartPtr<DoFDistribution> dd,
						  number refTol,
						  number coarsenTol,
						  int minLevel,
						  int maxLevel,
						  bool refTopLvlOnly = false)
{
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	int numMarkedRefine = 0, numMarkedCoarse = 0;
	const_iterator iter = dd->template begin<TElem>();
	const_iterator iterEnd = dd->template end<TElem>();
	const MultiGrid* mg = dd->multi_grid().get();
	int topLvl = 0;
	if(mg)
		topLvl = (int)mg->top_level();
	else
		refTopLvlOnly = false;

//	loop elements for marking
	for(; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

	//	marks for refinement
		if((refTol >= 0)
			&& (aaError[elem] > refTol)
			&& (dd->multi_grid()->get_level(elem) < maxLevel)
			&& ((!refTopLvlOnly) || (mg->get_level(elem) == topLvl)))
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
	if(pcl::NumProcs() > 1)
	{
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

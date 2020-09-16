/*
 * Copyright (c) 2015:  G-CSC, Goethe University Frankfurt
 * Author: Arne Nägel
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


#ifndef __H__UG_DISC__ERROR_ELEM_MARKING_STRATEGY__
#define __H__UG_DISC__ERROR_ELEM_MARKING_STRATEGY__

#include "lib_grid/algorithms/debug_util.h"  // for ElementDebugInfo
#include "lib_grid/multi_grid.h"
#include "lib_grid/refinement/refiner_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "error_indicator_util.h"

namespace ug{

/// Abstract base class for element marking (in adaptive refinement)
template <typename TDomain> class IElementMarkingStrategy;


/// This class encapsulates the multi-grid attachments for error estimation
/** Purpose: replaces direct access to 'm_aaError' etc. */
template <typename TDomain>
class IMultigridElementIndicators
{
public:
	///	world dimension
	static const int dim = TDomain::dim;

	typedef typename domain_traits<dim>::element_type elem_type;
	typedef Attachment<number> error_attachment_type;
	typedef MultiGrid::AttachmentAccessor<elem_type, error_attachment_type > attachment_accessor_type;

	/// CTOR
	IMultigridElementIndicators() {}

	/// DTOR
	~IMultigridElementIndicators()
	{
		detach_indicators();
	}

	/// Attach error indicator to multigrid
	void attach_indicators(SmartPtr<MultiGrid> pMG)
	{
		typedef typename domain_traits<dim>::element_type elem_type;

		if (!m_pMG->has_attachment<elem_type>(m_aError))
			pMG->template attach_to_dv<elem_type>(m_aError, -1.0);  // attach with default value
		m_pMG = pMG;
		m_aaError = attachment_accessor_type(*m_pMG, m_aError);
	}

	/// Detach error indicator from multigrid
	void detach_indicators()
	{
		if (m_pMG.invalid()) return; // no elements attached
		if (m_pMG->has_attachment<elem_type>(m_aError))
			m_pMG->template detach_from<elem_type>(m_aError);
	}


	/// returns error indicator value
	number& error(typename attachment_accessor_type::atraits::ConstElemPtr pElem)
	{ return this->m_aaError[pElem]; }

	/// returns error indicator value
	const number& error(typename attachment_accessor_type::atraits::ConstElemPtr pElem) const
	{ return this->m_aaError[pElem]; }


	friend class IElementMarkingStrategy<TDomain>;


	/// TODO: remove this function
	/// (mbreit: no, please leave it, it is very useful, at least with const access)
	attachment_accessor_type& errors()
	{ return m_aaError; }

protected:
	SmartPtr<MultiGrid> m_pMG; 		 // multigrid
	error_attachment_type m_aError;  // holding 'number' attachments
	attachment_accessor_type m_aaError;
};




/**
 *
 */
template <typename TDomain>
class IElementMarkingStrategy
{
public:
	///	world dimension
	static const int dim = TDomain::dim;

	/// element type to be marked
	typedef typename domain_traits<dim>::element_type elem_type;
	typedef typename Grid::AttachmentAccessor<elem_type, ug::Attachment<number> > elem_accessor_type;

	IElementMarkingStrategy()
	: m_latest_error(-1), m_latest_error_per_elem_max(-1), m_latest_error_per_elem_min(-1)
	{}
	virtual ~IElementMarkingStrategy(){};


	/// This function marks all elements
	void mark(IMultigridElementIndicators<TDomain>& mgElemIndicators,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
	{
		// forward call
		mark(mgElemIndicators.errors(), refiner, dd);
	};

	protected:
		/// DEPRECATED:
		virtual void mark(elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd) = 0;

		number m_latest_error;
		number m_latest_error_per_elem_max;
		number m_latest_error_per_elem_min;

	public:
		number global_estimated_error() const {return m_latest_error;}
		number global_estimated_error_per_elem_max() const {return m_latest_error_per_elem_max;}
		number global_estimated_error_per_elem_min() const {return m_latest_error_per_elem_min;}
};

/// M. Breit's standard refinement strategy
template <typename TDomain>
class StdRefinementMarkingStrategy : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	StdRefinementMarkingStrategy(number tol, int max_level)
	: m_tol(tol), m_max_level(max_level) {};

	void set_tolerance(number tol) {m_tol = tol;}
	void set_max_level(int max_level) {m_max_level = max_level;}

protected:
	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

	number m_tol;
	int m_max_level;
};

template <typename TDomain>
void StdRefinementMarkingStrategy<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	MarkElementsForRefinement<typename base_type::elem_type>(aaError, refiner, dd, m_tol, m_max_level);
}



/// mark everything if error too high and refinement allowed
template <typename TDomain>
class GlobalMarking : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	GlobalMarking(number tol, size_t max_level)
	: m_tol(tol), m_max_level(max_level) {};

	void set_tolerance(number tol) {m_tol = tol;}
	void set_max_level(size_t max_level) {m_max_level = max_level;}

protected:
	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

protected:

	number m_tol;
	size_t m_max_level;
};

template <typename TDomain>
void GlobalMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	number minElemErr;
	number maxElemErr;
	number errTotal;
	size_t numElem;

	ComputeMinMaxTotal(aaError, dd, minElemErr, maxElemErr, errTotal, numElem);
	if (errTotal <= m_tol || dd->multi_grid()->num_levels() > m_max_level)
		return;

	typedef typename DoFDistribution::traits<typename base_type::elem_type>::const_iterator const_iterator;

	const_iterator iter = dd->template begin<typename base_type::elem_type>();
	const_iterator iterEnd = dd->template end<typename base_type::elem_type>();

//	loop elements for marking
	for (; iter != iterEnd; ++iter)
		refiner.mark(*iter, RM_REFINE);
}



/// M. Breit's standard coarsening strategy
template <typename TDomain>
class StdCoarseningMarkingStrategy : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	StdCoarseningMarkingStrategy(number tol)
		: m_tol(tol), m_safety(8.0), m_minLvl(0) {}

	StdCoarseningMarkingStrategy(number tol, number safety)
		: m_tol(tol), m_safety(safety), m_minLvl(0) {}

	StdCoarseningMarkingStrategy(number tol, number safety, int minLvl)
		: m_tol(tol), m_safety(safety), m_minLvl(minLvl) {}

	void set_tolerance(number tol) {m_tol = tol;}
	void set_safety_factor(number safety) {m_safety = safety;}

protected:
	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

protected:
	number m_tol;
	number m_safety;
	int m_minLvl;
};

template <typename TDomain>
void StdCoarseningMarkingStrategy<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	MarkElementsForCoarsening<typename base_type::elem_type>(aaError, refiner, dd, m_tol, m_safety, m_minLvl);
}


/* Generate an ordered list of \eta_k^2 (in descending order, ie. largest first) */
template<class TElem>
number CreateListOfElemWeights(
		Grid::AttachmentAccessor<TElem, ug::Attachment<number> > &aaError,
		typename DoFDistribution::traits<TElem>::const_iterator iterBegin,
		const typename DoFDistribution::traits<TElem>::const_iterator iterEnd,
		std::vector<double> &etaSq)
{
	number localErrSq=0;
	typename DoFDistribution::traits<TElem>::const_iterator iter;
	size_t i=0;
	for (iter = iterBegin; iter != iterEnd; ++iter)
	{
		const double elemErrSq = aaError[*iter];

		//	If no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (elemErrSq < 0) continue;

		etaSq[i++]  = elemErrSq;
		localErrSq += elemErrSq;
	}

	// Sort in descending order (using default comparison).
	std::sort (etaSq.begin(), etaSq.end(), std::greater<double>());
	return localErrSq;
};



/* Generate an ordered list of \eta_k^2 (in descending order, ie. largest first) */
template<class TElem>
number CreateSortedListOfElems(
		Grid::AttachmentAccessor<TElem, ug::Attachment<number> > &aaError,
		typename DoFDistribution::traits<TElem>::const_iterator iterBegin,
		const typename DoFDistribution::traits<TElem>::const_iterator iterEnd,
		std::vector< std::pair<double, TElem*> > &etaSqList)
{

	//typedef typename std::pair<double, TElem*> TPair;
	typename DoFDistribution::traits<TElem>::const_iterator iter;

	number localErrSq=0;
	size_t i=0;
	for (iter = iterBegin; iter != iterEnd; ++iter)
	{
		const double elemErrSq = aaError[*iter];

		//	If no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (elemErrSq < 0) continue;

		etaSqList[i++]  = std::make_pair<>(elemErrSq, *iter);
		localErrSq += elemErrSq;
	}

	// Sort in descending order (using default comparison).
	std::sort (etaSqList.begin(), etaSqList.end());
	return localErrSq;
};



/// Refine as many largest-error elements as necessary to reach tolerance.
/// The expected tolerance is calculated using a user-given expected error reduction factor.
template <typename TDomain>
class ExpectedErrorMarkingStrategy : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	ExpectedErrorMarkingStrategy(number tol, int max_level, number safetyFactor, number expectedReductionFactor)
	: m_tol(tol), m_max_level(max_level), m_safety(safetyFactor), m_expRedFac(expectedReductionFactor) {};

	void set_tolerance(number tol) {m_tol = tol;}
	void set_max_level(int max_level) {m_max_level = max_level;}
	void set_safety_factor(number safetyFactor) {m_safety = safetyFactor;}
	void set_expected_reduction_factor(number expectedReductionFactor) {m_expRedFac = expectedReductionFactor;}

	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

protected:
	void merge_sorted_lists
	(
		std::vector<std::pair<number, int> >::iterator beginFirst,
		std::vector<std::pair<number, int> >::iterator beginSecond,
		std::vector<std::pair<number, int> >::iterator end
	);

protected:
	number m_tol;
	int m_max_level;
	number m_safety;
	number m_expRedFac;
};


template <typename TDomain>
struct ElemErrorSortDesc
{
	typedef typename domain_traits<TDomain::dim>::element_type elem_type;
	typedef typename Grid::AttachmentAccessor<elem_type, ug::Attachment<number> > error_accessor_type;
	ElemErrorSortDesc(const error_accessor_type& aaErr)
	: m_aaErr(aaErr) {}

	bool operator()(const elem_type* elem1, const elem_type* elem2)
	{
		return m_aaErr[elem1] > m_aaErr[elem2];
	}

	const error_accessor_type& m_aaErr;
};


template <typename TDomain>
void ExpectedErrorMarkingStrategy<TDomain>::merge_sorted_lists
(
	std::vector<std::pair<number, int> >::iterator beginFirst,
	std::vector<std::pair<number, int> >::iterator beginSecond,
	std::vector<std::pair<number, int> >::iterator end
)
{
	const size_t nVal = std::distance(beginFirst, end);
	std::vector<std::pair<number, int> > sorted;
	sorted.reserve(nVal);

	std::vector<std::pair<number, int> >::iterator it1 = beginFirst;
	std::vector<std::pair<number, int> >::iterator it2 = beginSecond;
	while (it1 != beginSecond && it2 != end)
	{
		if (it1->first > it2->first)
		{
			sorted.push_back(*it1);
			++it1;
		}
		else
		{
			sorted.push_back(*it2);
			++it2;
		}
	}

	for (; it1 != beginSecond; ++it1)
		sorted.push_back(*it1);
	for (; it2 != end; ++it2)
		sorted.push_back(*it2);

	size_t i = 0;
	for (it1 = beginFirst; it1 != end; ++it1, ++i)
		*it1 = sorted[i];
}



template <typename TDomain>
void ExpectedErrorMarkingStrategy<TDomain>::mark
(
	typename base_type::elem_accessor_type& aaErrorSq,
	IRefiner& refiner,
	ConstSmartPtr<DoFDistribution> dd
)
{
	typedef typename base_type::elem_type TElem;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	// create vector of local element (squared) errors
	const_iterator iter = dd->template begin<TElem>();
	const const_iterator iterEnd = dd->template end<TElem>();
	std::vector<TElem*> elemVec;
	number locError = 0.0;
	for (; iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		if (aaErrorSq[elem] < 0)
			continue;

		if (aaErrorSq[elem] != aaErrorSq[elem])
		{
			UG_LOG_ALL_PROCS("NaN error indicator on " << ElementDebugInfo(*refiner.grid(), elem));
			continue;
		}

		locError += aaErrorSq[elem];

		if (dd->multi_grid()->get_level(elem) < m_max_level)
			elemVec.push_back(elem);
	}
	const size_t nLocalElem = elemVec.size();

	// sort vector of elements locally (descending)
	ElemErrorSortDesc<TDomain> eeSort(aaErrorSq);
	std::sort(elemVec.begin(), elemVec.end(), eeSort);

	// create sorted vector of errors
	std::vector<number> etaSq(nLocalElem);
	for (size_t i = 0; i < nLocalElem; ++i)
		etaSq[i] = aaErrorSq[elemVec[i]];

#ifdef UG_PARALLEL
	const int nProcs = pcl::NumProcs();
	if (nProcs > 1)
	{
		// calculate global error
		pcl::ProcessCommunicator pc;
		const int rootProc = pcl::NumProcs() - 1;
		const number globError = pc.allreduce(locError, PCL_RO_SUM);
		UG_LOGN("  +++ Element errors: sumEtaSq = " << globError << ".");

		// construct global sorted vector of elem errors
		// gather on root proc
		// choose root as the highest-rank proc (proc 0 has least free memory)
		std::vector<number> rcvBuf;
		std::vector<int> vSizes;
		pc.gatherv(rcvBuf, etaSq, rootProc, &vSizes);

		// merge of pre-sorted local vectors
		std::vector<int> nRefineElems;
		size_t globNumRefineElems = 0;
		if (pcl::ProcRank() == rootProc)
		{
			// associate each error value with the proc it belongs to
			// and calculate offsets at the same time
			std::vector<size_t> vOffsets(nProcs, 0);
			const size_t nGlobElem = rcvBuf.size();
			std::vector<std::pair<number, int> > vGlobErrors(nGlobElem);
			size_t e = 0;
			for (int p = 0; p < nProcs; ++p)
			{
				const size_t szp = vSizes[p];
				for (size_t i = 0; i < szp; ++i, ++e)
				{
					vGlobErrors[e].first = rcvBuf[e];
					vGlobErrors[e].second = p;
				}
				if (p < nProcs-1)
					vOffsets[p+1] = e;
			}

			// free rcv buffer
			{
				std::vector<number> tmp;
				tmp.swap(rcvBuf);
			}

			// merge pre-sorted proc error lists ((log p) * N)
			size_t nLists = 2;
			while (true)
			{
				const size_t nMerge = nProcs / std::min(nLists, (size_t) nProcs);
				for (size_t m = 0; m < nMerge; ++m)
				{
					std::vector<std::pair<number, int> >::iterator beginFirst =
						vGlobErrors.begin() + vOffsets[m*nLists];
					std::vector<std::pair<number, int> >::iterator beginSecond =
						vGlobErrors.begin() + vOffsets[m*nLists + nLists/2];
					std::vector<std::pair<number, int> >::iterator endSecond =
						(m+1)*nLists < (size_t) nProcs ? vGlobErrors.begin() + vOffsets[(m+1)*nLists] : vGlobErrors.end();

					merge_sorted_lists(beginFirst, beginSecond, endSecond);
				}

				if (nLists >= (size_t) nProcs)
					break;

				nLists = nLists << 1;
			}

			// decide how many elements on each proc have to be refined
			const number requiredReduction = globError > m_tol ? globError - m_safety*m_tol : 0.0;
			number red = 0.0;
			nRefineElems.resize(nProcs, 0);
			for (; red < requiredReduction && globNumRefineElems < nGlobElem; ++globNumRefineElems)
			{
				red += (1.0-m_expRedFac) * vGlobErrors[globNumRefineElems].first;
				++nRefineElems[vGlobErrors[globNumRefineElems].second];
			}
		}

		// tell each proc how many of their elements are to be marked
		int nRefElems = 0;
		pc.scatter(GetDataPtr(nRefineElems), 1, PCL_DT_INT, &nRefElems, 1, PCL_DT_INT, rootProc);

		pc.broadcast(globNumRefineElems, rootProc);

		// mark for refinement
		UG_COND_THROW((size_t) nRefElems > nLocalElem, "More elements supposedly need refinement here ("
			<< nRefElems << ") than are refineable (" << nLocalElem << ").");
		for (size_t i = 0; i < (size_t) nRefElems; ++i)
			refiner.mark(elemVec[i], RM_REFINE);

		if (globNumRefineElems)
		{
			UG_LOGN("  +++ Marked for refinement: " << globNumRefineElems << " elements");
		}
		else
		{
			UG_LOGN("  +++ No refinement necessary.");
		}
	}
	else
	{
#endif

	// plan refinement of elements until expected error reduction is enough
	// to get below tolerance (or all elements are marked)
	const number requiredReduction = locError - m_safety*m_tol;
	UG_LOGN("  +++ Element errors: sumEtaSq = " << locError << ".");
	number red = 0.0;
	size_t i = 0;
	for (; red < requiredReduction && i < nLocalElem; ++i)
	{
		red += m_expRedFac * etaSq[i];
		refiner.mark(elemVec[i], RM_REFINE);
	}
	if (i)
	{
		UG_LOGN("  +++ Marked for refinement: " << i << " elements.");
	}
	else
	{
		UG_LOGN("  +++ No refinement necessary.");
	}

#ifdef UG_PARALLEL
	}
#endif
}



/// Marks elements with \eta_K >= \theta  \max_{K'} \eta_{K'} for refinement
//! (cf. Verfuerth script)
template <typename TDomain>
class MaximumMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	MaximumMarking(number theta=1.0)
	: m_theta(theta), m_theta_min(0.0), m_eps(0.01), m_max_level(100), m_min_level(0) {};
	MaximumMarking(number theta, number eps)
	: m_theta(theta), m_theta_min(0.0), m_eps (eps), m_max_level(100), m_min_level(0) {};
	MaximumMarking(number theta_max, number theta_min, number eps)
	: m_theta(theta_max), m_theta_min(theta_min), m_eps (eps), m_max_level(100), m_min_level(0) {};

protected:
	void mark(typename base_type::elem_accessor_type& aaErrorSq, IRefiner& refiner, ConstSmartPtr<DoFDistribution> dd);

public:
	void set_max_level(int lvl) {m_max_level = lvl;}
	void set_min_level(int lvl) {m_min_level = lvl;}

protected:

	number m_theta, m_theta_min;
	number m_eps;
	int m_max_level, m_min_level;
};



template <typename TDomain>
void MaximumMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaErrorSq,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename base_type::elem_type TElem;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	// compute minimal/maximal/ total error and number of elements

	number minElemErr, minElemErrLocal;
	number maxElemErr, maxElemErrLocal;
	number errTotal, errLocal;
	size_t numElem, numElemLocal;

	ComputeMinMax(aaErrorSq, dd, minElemErr, maxElemErr, errTotal, numElem,
				minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);

	this->m_latest_error = sqrt(errTotal);
	this->m_latest_error_per_elem_max = maxElemErr;
	this->m_latest_error_per_elem_min = minElemErr;

	// init iterators
	const_iterator iter;
	const const_iterator iterEnd = dd->template end<TElem>();

	// determine (local) number of excess elements
	const size_t ndiscard = (size_t) (numElemLocal*m_eps); // TODO: on every process?
	UG_LOG("  +++ MaximumMarking: Found max "<<  maxElemErr << ", ndiscard="<<ndiscard<<".\n");

	if (numElemLocal > 0)
	{
		// create sorted array of elem weights
		std::vector<double> etaSq;
		etaSq.resize(numElemLocal);
		CreateListOfElemWeights<TElem>(aaErrorSq,dd->template begin<TElem>(), iterEnd, etaSq);
		UG_ASSERT(numElemLocal==etaSq.size(), "Huhh: number of elements does not match!");
		UG_ASSERT(numElemLocal > ndiscard, "Huhh: number of elements does not match!");
		maxElemErr = etaSq[ndiscard];
	}

	// compute parallel threshold
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number maxElemErrLocal = maxElemErr;
		maxElemErr = com.allreduce(maxElemErrLocal, PCL_RO_MAX);
		UG_LOG("  +++ Maximum for refinement: " << maxElemErr << ".\n");
	}
#else
	UG_LOG("  +++ Skipping " << ndiscard << " elements; new (parallel) max." << maxElemErr << ".\n");
#endif

	// refine all element above threshold
	const number top_threshold = maxElemErr*m_theta;
	const number bot_threshold = maxElemErr*m_theta_min;
	UG_LOG("  +++ Refining elements, if error > " << maxElemErr << "*" << m_theta <<
			" = " << top_threshold << ".\n");
	UG_LOG("  +++ Coarsening elements, if error < " << maxElemErr << "*" << m_theta_min <<
				" = " << bot_threshold << ".\n");

	//	reset counter
	std::size_t numMarkedRefine = 0;
	std::size_t numMarkedCoarse = 0;

	//	loop elements for marking
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		//	if no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (aaErrorSq[elem] < 0) continue;

		//	mark for refinement
		if ((aaErrorSq[elem] > top_threshold) && (dd->multi_grid()->get_level(elem) <= m_max_level))
		{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
		}

		//	mark for coarsening
		if ((aaErrorSq[elem] < bot_threshold) && (dd->multi_grid()->get_level(elem) >= m_min_level))
		{
				refiner.mark(elem, RM_COARSEN);
				numMarkedCoarse++;
		}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::size_t numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		UG_LOG("  +++ MaximumMarking: Marked for refinement: " << numMarkedRefine << " ("<< numMarkedRefineLocal << ") elements.\n");
	}
	else
#endif
	{ UG_LOG("  +++ MaximumMarking: refinement: " << numMarkedRefine << ", coarsening" << numMarkedCoarse <<  " elements.\n") }


}



/// Marks element with smallest \eta_i^2
//! (cf. Verfuerth script)
template <typename TDomain>
class APosterioriCoarsening : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	APosterioriCoarsening(number theta=0.1)
	: m_theta(theta), m_max_level(100), m_min_level(0) {};
protected:
	void mark(typename base_type::elem_accessor_type& aaErrorSq, IRefiner& refiner, ConstSmartPtr<DoFDistribution> dd);
public:
	void set_max_level(int lvl) {m_max_level = lvl;}
	void set_min_level(int lvl) {m_min_level = lvl;}

protected:

	number m_theta;
	int m_max_level, m_min_level;
};



template <typename TDomain>
void APosterioriCoarsening<TDomain>::mark(typename base_type::elem_accessor_type& aaErrorSq,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename base_type::elem_type TElem;
	//typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;
	typedef typename std::pair<double, TElem*> TPair;
	typedef typename std::vector<TPair> TPairVector;

	// Compute minimal/maximal/ total error and number of elements.
	number minElemErrSq, minElemErrSqLocal;
	number maxElemErrSq, maxElemErrSqLocal;
	number errSqTotal, errSqLocal;
	size_t numElem, numElemLocal;

	size_t numCoarsened=0;


	ComputeMinMax(aaErrorSq, dd, minElemErrSq, maxElemErrSq, errSqTotal, numElem,
				minElemErrSqLocal, maxElemErrSqLocal, errSqLocal, numElemLocal);

	this->m_latest_error = sqrt(errSqTotal);
	this->m_latest_error_per_elem_max = maxElemErrSq;
	this->m_latest_error_per_elem_min = minElemErrSq;



	if (numElemLocal > 0)
	{
		// Create sorted array of elem weights.
		TPairVector etaSqVec;
		etaSqVec.resize(numElemLocal);
		CreateSortedListOfElems<TElem>(aaErrorSq,dd->template begin<TElem>(),  dd->template end<TElem>(), etaSqVec);


		UG_ASSERT(numElemLocal==etaSqVec.size(), "Huhh: number of elements does not match!");

		{

			const double mySumSq = errSqLocal*m_theta;
			double localSumSq = 0.0;

			typename TPairVector::const_iterator iterEnd = etaSqVec.end();

			for (typename TPairVector::iterator iter = etaSqVec.begin();
				(iter != iterEnd) && (localSumSq < mySumSq); ++iter)
			{

				localSumSq += iter->first;
				if (localSumSq < mySumSq) {
					refiner.mark(iter->second, RM_COARSEN);
					numCoarsened++;
				}

			}
			UG_LOG("  +++ APosterioriCoarsening: localSumSq = "<< localSumSq << " >" << mySumSq << std::endl);

		}



	}


#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::size_t numCoarsenedLocal = numCoarsened;
		numCoarsened = com.allreduce(numCoarsenedLocal, PCL_RO_SUM);
		UG_LOG("  +++ APosterioriCoarsening: Marked for coarsening: " << numCoarsened << " ("<< numCoarsenedLocal << ") elements.\n");
	}
	else
#endif
	{ UG_LOG("  +++ APosterioriCoarsening: coarsening" << numCoarsened <<  " elements.\n") }


}


//! marks elements above a certain fraction of the maximum
/*! Cf. Verfuerth' scriptum */
template <typename TDomain>
class EquilibrationMarkingStrategy : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	EquilibrationMarkingStrategy(number theta_top=0.9)
	: m_theta_top(theta_top), m_theta_bot(0.0) {} //, m_max_level(100) {};
	EquilibrationMarkingStrategy(number theta_top, number theta_bot)
	: m_theta_top(theta_top), m_theta_bot(theta_bot) {} ;// , m_max_level(100) {};


protected:
	void mark(typename base_type::elem_accessor_type& aaErrorSq,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);

	number m_theta_top;  			// refine at least a certain fraction,  0.0 <= m_theta_top <= 1.0
	number m_theta_bot;  			// refine at least a certain fraction,  0.0 <= m_theta_bot <= 1.0
	// int m_max_level;
};





template <typename TDomain>
void EquilibrationMarkingStrategy<TDomain>::mark(typename base_type::elem_accessor_type& aaErrorSq,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename base_type::elem_type TElem;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	// compute minimal/maximal/total error and number of elements
	number minElemErr, minElemErrLocal;
	number maxElemErr, maxElemErrLocal;
	number errTotal= 0.0, errLocal= 0.0;
	size_t numElem= 0, numElemLocal;

	// compute weights
	ComputeMinMax(aaErrorSq, dd, minElemErr, maxElemErr, errTotal, numElem,
					minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);

	// init iterators
	const const_iterator iterEnd = dd->template end<TElem>();
	const_iterator iter;

	// create and fill sorted array of $\etaSq^2_i$ for all (local) elements
	std::vector<double> etaSq;
	etaSq.resize(numElemLocal);
	CreateListOfElemWeights<TElem>(aaErrorSq,dd->template begin<TElem>(), iterEnd, etaSq);
	UG_ASSERT(numElemLocal==etaSq.size(), "Huhh: number of elements does not match!");

	// compute thresholds
	UG_ASSERT( ((m_theta_top>=0.0) && (m_theta_top<=1.0)), "Huhh: m_theta_top invalid!");
	UG_ASSERT( ((m_theta_bot>=0.0) && (m_theta_bot<=1.0)), "Huhh: m_theta_top invalid!");
	UG_ASSERT( (m_theta_top>m_theta_bot), "Huhh: m_theta_top invalid!");

	// discard a fraction of elements
	// a) largest elements
	typename std::vector<double>::const_iterator top = etaSq.begin();
	for (number sumSq=0.0;
		(sumSq<m_theta_top*errLocal) && (top !=(etaSq.end()-1));
		++top) { sumSq += *top; }
	number top_threshold = (top != etaSq.begin()) ? (*top + *(top-1))/2.0 : *top;

	// a) smallest elements
	typename std::vector<double>::const_iterator bot = etaSq.end()-1;
		for (number sumSq=0.0;
			(sumSq<m_theta_bot*errLocal) && (bot !=etaSq.begin() );
			--bot) { sumSq += *bot; }
	number bot_threshold = (bot != (etaSq.end()-1)) ? (*bot + *(bot+1))/2.0 : 0.0;

	UG_LOG("  +++  error = "<<  errLocal << std::endl);
	UG_LOG("  +++  top_threshold= "<< top_threshold <<"( "<< top - etaSq.begin() << " cells)" << std::endl);
	UG_LOG("  +++  bot_threshold= "<< bot_threshold <<"( "<< etaSq.end()-bot << " cells)" << std::endl);

	//	mark elements with maximal contribution
	size_t numMarkedRefine = 0;
	size_t numMarkedCoarse = 0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{

		const double elemErr = aaErrorSq[*iter];		// get element
		if (elemErr < 0) continue;					// skip invalid

		if (elemErr > top_threshold)
		{
			refiner.mark(*iter, RM_REFINE);
			numMarkedRefine++;
		}

		if (elemErr < bot_threshold)
		{
			refiner.mark(*iter, RM_COARSEN);
			numMarkedCoarse++;
		}
	}



	UG_LOG("  +++ EquilibrationMarkingStrategy: Marked for refinement: "<<  numMarkedRefine << ", " << numMarkedCoarse << " elements.\n");

}


/// Marks elements above \f$ \theta * (\mu + width * \sigma) \f$ for refinement
//! where \f$ \mu = E[\eta^2], \sigma^2 = Var[\eta^2] \f$
template <typename TDomain>
class VarianceMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	VarianceMarking(number theta) : m_theta(theta), m_width(3.0), m_max_level(100) {};
	VarianceMarking(number theta, number width) : m_theta(theta), m_width (width), m_max_level(100) {};


protected:
	void mark(typename base_type::elem_accessor_type& aaError2,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);

	number m_theta;
	number m_width;
	int m_max_level;
};

template <typename TDomain>
void VarianceMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaError2,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename base_type::elem_type TElem;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	// compute minimal/maximal/ total error and number of elements

	number minElemErr, minElemErrLocal;
	number maxElemErr, maxElemErrLocal;
	number errTotal, errLocal;
	size_t numElem, numElemLocal;

	ComputeMinMax(aaError2, dd, minElemErr, maxElemErr, errTotal, numElem,
				minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);

	this->m_latest_error = sqrt(errTotal);
	this->m_latest_error_per_elem_max = maxElemErr;
	this->m_latest_error_per_elem_min = minElemErr;

//	number elemMean =  sqrt(errTotal) / numElem;
	number elemMean =  errTotal / numElem;

	UG_LOG("  +++ VarianceMarking: Mean error : " << elemMean << " on "<< numElem << " elements.\n");

	// init iterators
	const_iterator iter;
	const const_iterator iterEnd = dd->template end<TElem>();

	number elemVar = 0.0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;
		number elemError = aaError2[elem];  // eta_i^2

		if (elemError < 0) continue;
		elemVar += (elemMean-elemError) * (elemMean-elemError);


	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number elemVarLocal = elemVar;
		elemVar = com.allreduce(elemVarLocal, PCL_RO_SUM);

	}
#endif
	UG_LOG("  +++ VarianceMarking: Est. variance (1) : " << elemVar << " on "<< numElem << " elements.\n");


	elemVar /= (numElem-1.0);
	UG_LOG("  +++ VarianceMarking: Est. variance (2): " << elemVar << " on "<< numElem << " elements.\n");


		// refine all element above threshold
	const number sigma = sqrt(elemVar);
	const number maxError = elemMean + sigma*m_width;
	UG_LOG("  +++ Refining elements if error greater " << sigma << "*" << m_width <<" + "<< elemMean <<
			" = " << maxError << ".\n");



	const number minErrToRefine = maxError*m_theta;
	UG_LOG("  +++ Refining elements if error greater " << maxError << "*" << m_theta <<
				" = " << minErrToRefine << ".\n");

	//	reset counter
	std::size_t numMarkedRefine = 0;

	//	loop elements for marking
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		//	if no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (aaError2[elem] < 0) continue;

		//	marks for refinement
		if (aaError2[elem] >= minErrToRefine)
			if (dd->multi_grid()->get_level(elem) <= m_max_level)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::size_t numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		UG_LOG("  +++ MaximumMarking: Marked for refinement: " << numMarkedRefine << " ("<< numMarkedRefineLocal << ") elements.\n");
	}
#else
	UG_LOG("  +++ MaximumMarking: Marked for refinement: " << numMarkedRefine << " elements.\n");
#endif


}


/// Marks elements above a certain threshold for refinement
/*! Threshold given by
 * \f$ \theta * (\mu + width * \sigma) \f$
 *  where
 *  \f$ \mu = E[\eta], \sigma^2 = Var[\eta] \f$*
 */
template <typename TDomain>
class VarianceMarkingEta : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	VarianceMarkingEta(number theta) :
		m_theta(theta), m_width(3.0), m_max_level(100), m_theta_coarse(0.0), m_min_level(0)
	{};
	VarianceMarkingEta(number theta, number width) :
		m_theta(theta), m_width (width), m_max_level(100), m_theta_coarse(0.0), m_min_level(0)
	{};
	VarianceMarkingEta(number theta, number width, number theta_coarse) :
			m_theta(theta), m_width (width), m_max_level(100), m_theta_coarse(theta_coarse), m_min_level(0)
	{};

	void init_refinement(number theta, int max_level)
	{m_theta = theta; m_max_level=max_level;}

	void init_coarsening(number theta, int min_level)
	{m_theta_coarse = theta; m_min_level=min_level;}

protected:
	void mark(typename base_type::elem_accessor_type& aaError2,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_theta;
	number m_width;
	int m_max_level;

	number m_theta_coarse;
	int m_min_level;
};

template <typename TDomain>
void VarianceMarkingEta<TDomain>::mark(typename base_type::elem_accessor_type& aaError2,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename base_type::elem_type TElem;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	// compute minimal/maximal/ total error and number of elements

	number minElemErr, minElemErrLocal;
	number maxElemErr, maxElemErrLocal;
	number errSum, errTotalSq, errLocal;
	size_t numElem, numElemLocal;

	const number elemMean = ComputeAvg(aaError2, dd, minElemErr, maxElemErr, errSum, errTotalSq, numElem,
									minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);


	this->m_latest_error = sqrt(errTotalSq);
	this->m_latest_error_per_elem_max = maxElemErr;
	this->m_latest_error_per_elem_min = minElemErr;

	UG_LOG("  +++ VarianceMarkingEta: error : "<< this->m_latest_error << " (meanEta : " << elemMean << " on "<< numElem << " elements).\n");

	// init iterators
	const_iterator iter;
	const const_iterator iterEnd = dd->template end<TElem>();

	number elemVar = 0.0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		if (aaError2[elem] < 0) continue;

		number elemError =  sqrt(aaError2[elem]); // eta_i

		elemVar += (elemMean-elemError) * (elemMean-elemError);
	}

	UG_LOG("  +++ VarianceMarkingEta: Est. variance (1) : " << (elemVar / (numElem -1.0)) << " on "<< numElem << " elements.\n");
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number elemVarLocal = elemVar;
		elemVar = com.allreduce(elemVarLocal, PCL_RO_SUM);

	}
#endif
	elemVar /= (numElem-1.0);
	UG_LOG("  +++ VarianceMarkingEta: Est. variance (2): " << elemVar << " on "<< numElem << " elements.\n");


	// refine all element above threshold
	const number sigma = sqrt(elemVar);
	const number maxError = elemMean + sigma*m_width;
	UG_LOG("  +++ Refining elements if error > " << elemMean << " + "<< sigma << "*" << m_width <<
			" = " << maxError << ".\n");

	const number elemErrorRefine = maxError*m_theta;
	UG_LOG("  +++ Refining elements if error > " << maxError << "*" << m_theta <<
				" = " << elemErrorRefine << ".\n");

	// coarsen all element not contributing significantly
	//const number minError = std::max(elemMean /*- sigma*m_width*/, 0.0);
	//UG_LOG("  +++ Coarsening elements if error greater " << elemMean << " - "<< sigma << "*" << m_width <<
	//			" = " << minError << ".\n");

	const number elemErrorCoarsen = elemMean*m_theta_coarse;
	UG_LOG("  +++ Coarsening elements if error < " << elemMean << "*" << m_theta_coarse <<
					" = " << elemErrorCoarsen << ".\n");



	//	reset counters
	std::size_t numMarkedRefine  = 0;
	std::size_t numMarkedCoarsen = 0;

	//	loop elements for marking
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		//	if no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (aaError2[elem] < 0) continue;

		const number elemError =  sqrt(aaError2[elem]); // eta_i

		//	mark for refinement
		if ((elemError >= elemErrorRefine) && (dd->multi_grid()->get_level(elem) <= m_max_level))
		{
			refiner.mark(elem, RM_REFINE);
			numMarkedRefine++;
		}

		//	mark for coarsening
		if ((elemError < elemErrorCoarsen) && (dd->multi_grid()->get_level(elem) >= m_min_level))
		{
			refiner.mark(elem, RM_COARSEN);
			numMarkedCoarsen++;
		}

	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::size_t numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		UG_LOG("  +++ VarianceMarkingEta: Marked for refinement: " << numMarkedRefine << " ("<< numMarkedRefineLocal << ") elements.\n");
	} else
#endif
	{
		UG_LOG("  +++ VarianceMarkingEta: Marked for refinement: " << numMarkedRefine << " elements.\n");
		UG_LOG("  +++ VarianceMarkingEta: Marked for coarsening: " << numMarkedCoarsen << " elements.\n");
	}
}


/// marks elements above a certain fraction of the maximum
//!
template <typename TDomain>
class MeanValueMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	MeanValueMarking(number theta, number factor) : m_theta(theta), m_factor (factor) {};

protected:
	void mark(typename base_type::elem_accessor_type& aaError,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_theta;
	number m_factor;
};

template <typename TDomain>
void MeanValueMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	typedef typename base_type::elem_type TElem;
	typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

	// compute minimal/maximal/ total error and number of elements

	number minElemErr, minElemErrLocal;
	number maxElemErrSq, maxElemErrLocal;
	number errTotalSq, errLocal;
	size_t numElem, numElemLocal;

	ComputeMinMax(aaError, dd, minElemErr, maxElemErrSq, errTotalSq, numElem,
				minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);


	// init iterators
	number avgElemErr = errTotalSq / numElem;
	UG_LOG("  +++ Global max: "<<  maxElemErrSq << ", Global min: "<< minElemErr <<".\n");
	UG_LOG("  +++ MeanValueMarking: Mean value : " << avgElemErr << " (Global error sum: "<<errTotalSq<<" on "<< numElem <<" elements).\n");

	// refine all element above threshold
	const number minThetaErrToRefine = maxElemErrSq*m_theta;
	const number minFactorAvgErrToRefine = avgElemErr*m_factor;
	UG_LOG("  +++ MeanValueMarking: Min theta error : "<<minThetaErrToRefine<<" (Global Max Error: " << maxElemErrSq<< " for theta: "<<m_theta<<").\n");
	UG_LOG("  +++ MeanValueMarking: Min factor avg error : "<<minFactorAvgErrToRefine<<" (Global Avg Error: " << avgElemErr<< " for factor: "<<m_factor<<").\n");

	const number minErrToRefine = std::min(minThetaErrToRefine, minFactorAvgErrToRefine);
	UG_LOG("  +++ MeanValueMarking: Refining if error >= : "<<minErrToRefine<<".\n");


	//	reset counter
	std::size_t numMarkedRefine = 0;

	//	loop elements for marking
	const_iterator iter;
	const const_iterator iterEnd = dd->template end<TElem>();
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		//	get element
		TElem* elem = *iter;

		//	if no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (aaError[elem] < 0) continue;

		//	marks for refinement
		if (aaError[elem] >= minErrToRefine)
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}
	}

#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::size_t numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		UG_LOG("  +++ MeanValueMarking: Marked for refinement: " << numMarkedRefine << " ("<< numMarkedRefineLocal << ") elements.\n");
	}
#else
	UG_LOG("  +++ MeanValueMarking: Marked for refinement: " << numMarkedRefine << " elements.\n");
#endif


}






/// marks elements above an absolute threshold (based on S. Reiter's idea)
template <typename TDomain>
class AbsoluteMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	AbsoluteMarking(number eta) : m_eta(eta), m_max_level(100) {};

protected:
	void mark(typename base_type::elem_accessor_type& aaError,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_eta;
	int m_max_level;
};

template <typename TDomain>
void AbsoluteMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
		typedef typename base_type::elem_type TElem;
		typedef typename DoFDistribution::traits<TElem>::const_iterator const_iterator;

		//	loop elements for marking
		const const_iterator iterEnd = dd->template end<TElem>();

		number lowerBound = m_eta*m_eta;

		number errTotal=0.0;

		size_t numMarkedRefine = 0;
		for (const_iterator iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
		{
		//	get element error
			TElem* elem = *iter;
			const number elemError = aaError[elem];

		// skip newly added
			if (elemError < 0) continue;
			errTotal += elemError;

		//	mark for refinement
			if ((elemError >= lowerBound) && (dd->multi_grid()->get_level(elem) <= m_max_level))
			{
				refiner.mark(elem, RM_REFINE);
				numMarkedRefine++;
			}
		}



#ifdef UG_PARALLEL
		if (pcl::NumProcs() > 1)
		{
			pcl::ProcessCommunicator com;
			std::size_t numMarkedRefineLocal = numMarkedRefine;
			number errLocal = errTotal;
			numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
			errTotal = com.allreduce(errLocal, PCL_RO_SUM);
		}
#endif

		this->m_latest_error = sqrt(errTotal);

		UG_LOG("  +++ AbsoluteMarking: Error**2 = "<<errTotal <<"; marked "<< numMarkedRefine << " elements for refinement "<<std::endl);

}



/// Mark surface layer for coarsening.
template <typename TDomain, typename TAlgebra>
void MarkForCoarsenening_SurfaceLayer(const GridFunction<TDomain, TAlgebra> &u, IRefiner& refiner){
		typedef typename domain_traits<TDomain::dim>::element_type TElem;
		ConstSmartPtr<DoFDistribution> spDD=u.dof_distribution();
		refiner.mark(spDD->begin<TElem>(), spDD->end<TElem>(), RM_COARSEN);
		refiner.coarsen();
}


}//	end of namespace

#endif


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

#include "lib_grid/multi_grid.h"
#include "lib_grid/algorithms/refinement/refiner_interface.h"
#include "lib_disc/dof_manager/dof_distribution.h"
#include "error_indicator_util.h"

namespace ug{

/// Abstract base class for element marking (in adaptive refinement)
template <typename TDomain>
struct IElementMarkingStrategy
{
	///	world dimension
	static const int dim = TDomain::dim;

	/// element type to be marked
	typedef typename domain_traits<dim>::element_type elem_type;
	typedef typename MultiGrid::AttachmentAccessor<elem_type, ug::Attachment<number> > elem_accessor_type;

	IElementMarkingStrategy() :
		m_latest_error(-1),
		m_latest_error_per_elem_max(-1),
		m_latest_error_per_elem_min(-1)
	{}
	virtual ~IElementMarkingStrategy(){};

	virtual void mark(elem_accessor_type& aaError,
			IRefiner& refiner,
			ConstSmartPtr<DoFDistribution> dd) = 0;

	protected:
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

	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

protected:

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
		: m_tol(tol), m_safety(8.0) {};

	StdCoarseningMarkingStrategy(number tol, number safety)
		: m_tol(tol), m_safety(safety) {};

	void set_tolerance(number tol) {m_tol = tol;}
	void set_safety_factor(number safety) {m_safety = safety;}

	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

protected:

	number m_tol;
	number m_safety;
};

template <typename TDomain>
void StdCoarseningMarkingStrategy<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	MarkElementsForCoarsening<typename base_type::elem_type>(aaError, refiner, dd, m_tol, m_safety);
}


template<class TElem>
number CreateListOfElemWeights(
		MultiGrid::AttachmentAccessor<TElem, ug::Attachment<number> > &aaError,
		typename DoFDistribution::traits<TElem>::const_iterator iterBegin,
		const typename DoFDistribution::traits<TElem>::const_iterator iterEnd,
		std::vector<double> &eta)
{
	number localErr=0;
	typename DoFDistribution::traits<TElem>::const_iterator iter; // = iterBegin;
	size_t i=0;
	for (iter = iterBegin; iter != iterEnd; ++iter)
	{
		const double elemErr = aaError[*iter];

			//	if no error value exists: ignore (might be newly added by refinement);
			//	newly added elements are supposed to have a negative error estimator
			if (elemErr < 0) continue;

			eta[i++]  = elemErr;
			localErr += elemErr;
	}

	// sort descending using default comparison
	std::sort (eta.begin(), eta.end(), std::greater<double>());
	return localErr;
};

/// marks elements above a certain fraction of the maximum
//!
template <typename TDomain>
class MaximumMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	MaximumMarking(number theta=1.0) : m_theta(theta), m_eps(0.01), m_max_level(100) {};
	MaximumMarking(number theta, number eps) : m_theta(theta), m_eps (eps), m_max_level(100) {};

	void mark(typename base_type::elem_accessor_type& aaError,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_theta;
	number m_eps;
	int m_max_level;
};

template <typename TDomain>
void MaximumMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
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

	ComputeMinMax(aaError, dd, minElemErr, maxElemErr, errTotal, numElem,
				minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);


	// init iterators
	const_iterator iter;
	const const_iterator iterEnd = dd->template end<TElem>();

	// determine (local) number of excess elements
	const size_t ndiscard = (size_t) (numElemLocal*m_eps); // TODO: on every process?
	UG_LOG("  +++ Found max "<<  maxElemErr << " ndiscard="<<ndiscard<<".\n");

	// Verfuerths strategy for skipping excess
/*	while (ndiscard)
	{
		// find elements that should be discarded
		int discard = 0;
		double newmax = 0.0;
		for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
		{
			const double elemErr = aaError[*iter];

			if ( elemErr >= maxElemErr) { discard++; }
			else { newmax = std::max(newmax, elemErr); }
		}

		// quit, if enough elements have been found
		if ( discard <  ndiscard) maxElemErr = newmax;
		else break;

	}*/

	if (numElemLocal > 0)
	{
		// create sorted array of elem weights
		std::vector<double> eta;
		eta.resize(numElemLocal);
		CreateListOfElemWeights<TElem>(aaError,dd->template begin<TElem>(), iterEnd, eta);
		UG_ASSERT(numElemLocal==eta.size(), "Huhh: number of elements does not match!");
		UG_ASSERT(numElemLocal > ndiscard, "Huhh: number of elements does not match!");
		maxElemErr = eta[ndiscard];

	}

	UG_LOG("  +++ MaximumMarking: Maximum for refinement: " << maxElemErr << " elements.\n");

	// compute parallel threshold
#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		number maxElemErrLocal = maxElemErr;
		maxElemErr = com.allreduce(maxElemErrLocal, PCL_RO_SUM);
		UG_LOG("  +++ MaximumMarking: Maximum for refinement: " << maxElemErr << ".\n");
		maxElemErr = maxElemErr / pcl::NumProcs();
		UG_LOG("  +++ MaximumMarking: Average for refinement: " << maxElemErr << ".\n");
	}
#else
	UG_LOG("  +++ Skipping " << ndiscard << " elements; new max." << maxElemErr << ".\n");
#endif

	// refine all element above threshold
	const number minErrToRefine = maxElemErr*m_theta;
	UG_LOG("  +++ Refining elements if error greater " << maxElemErr << "*" << m_theta <<
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
		if (aaError[elem] < 0) continue;

		//	marks for refinement
		if (aaError[elem] >= minErrToRefine)
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

/// marks elements above $\theta * (\mu + width * \sigma)$
//! where $\mu = E[\eta^2], \sigma^2 = Var[\eta^2]$
template <typename TDomain>
class VarianceMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	VarianceMarking(number theta) : m_theta(theta), m_width(3.0), m_max_level(100) {};
	VarianceMarking(number theta, number width) : m_theta(theta), m_width (width), m_max_level(100) {};

	void mark(typename base_type::elem_accessor_type& aaError,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_theta;
	number m_width;
	int m_max_level;
};

template <typename TDomain>
void VarianceMarking<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
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

	ComputeMinMax(aaError, dd, minElemErr, maxElemErr, errTotal, numElem,
				minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);

	this->m_latest_error = errTotal;
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
		number elemError = aaError[elem];  // eta_i^2

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
		if (aaError[elem] < 0) continue;

		//	marks for refinement
		if (aaError[elem] >= minErrToRefine)
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


/// marks elements above $\theta * (\mu + width * \sigma)$
//! where $\mu = E[\eta^2], \sigma^2 = Var[\eta^2]$
template <typename TDomain>
class VarianceMarkingEta : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	VarianceMarkingEta(number theta) : m_theta(theta), m_width(3.0), m_max_level(100) {};
	VarianceMarkingEta(number theta, number width) : m_theta(theta), m_width (width), m_max_level(100) {};

	void mark(typename base_type::elem_accessor_type& aaError,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_theta;
	number m_width;
	int m_max_level;
};

template <typename TDomain>
void VarianceMarkingEta<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
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

	const number elemMean = ComputeAvg(aaError, dd, minElemErr, maxElemErr, errSum, errTotalSq, numElem,
									minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);


	this->m_latest_error = sqrt(errTotalSq);
	this->m_latest_error_per_elem_max = maxElemErr;
	this->m_latest_error_per_elem_min = minElemErr;

//	number elemMean =  sqrt(errSum) / numElem;
	//number elemMean =  errSum / numElem;

	UG_LOG("  +++ VarianceMarkingEta: Mean error : " << elemMean << " on "<< numElem << " elements.\n");

	// init iterators
	const_iterator iter;
	const const_iterator iterEnd = dd->template end<TElem>();

	number elemVar = 0.0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		TElem* elem = *iter;

		if (aaError[elem] < 0) continue;

		number elemError =  sqrt(aaError[elem]); // eta_i

		elemVar += (elemMean-elemError) * (elemMean-elemError);


	}

	UG_LOG("  +++ VarianceMarkingEta: Est. variance (1) : " << elemVar << " on "<< numElem << " elements.\n");
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
		if (aaError[elem] < 0) continue;

		const number elemError =  sqrt(aaError[elem]); // eta_i

		//	marks for refinement
		if (elemError >= minErrToRefine)
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
		UG_LOG("  +++ VarianceMarkingEta: Marked for refinement: " << numMarkedRefine << " ("<< numMarkedRefineLocal << ") elements.\n");
	}
	else
		UG_LOG("  +++ VarianceMarkingEta: Marked for refinement: " << numMarkedRefine << " elements.\n");
#else
	UG_LOG("  +++ VarianceMarkingEta: Marked for refinement: " << numMarkedRefine << " elements.\n");
#endif


}


/// marks elements above a certain fraction of the maximum
//!
template <typename TDomain>
class MeanValueMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	MeanValueMarking(number theta, number factor) : m_theta(theta), m_factor (factor) {};

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




// marks elements above a certain fraction of the maximum
template <typename TDomain>
class EquilibrationMarkingStrategy : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	EquilibrationMarkingStrategy(number theta=1.0) : m_theta(theta), m_eps(0.01), m_max_level(100) {};
	EquilibrationMarkingStrategy(number theta, number eps) : m_theta(theta), m_eps (eps), m_max_level(100) {};

	void mark(typename base_type::elem_accessor_type& aaError,
					IRefiner& refiner,
					ConstSmartPtr<DoFDistribution> dd);
protected:

	number m_theta;
	number m_eps;  			// refine at least a certain fraction,  0.0 <= m_eps <= 1.0
	int m_max_level;
};





template <typename TDomain>
void EquilibrationMarkingStrategy<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
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

	// ComputeMinMaxTotal(aaError, dd, minElemErr, maxElemErr, errTotal, numTotalElem);
	ComputeMinMax(aaError, dd, minElemErr, maxElemErr, errTotal, numElem,
					minElemErrLocal, maxElemErrLocal, errLocal, numElemLocal);

	number localErr = 0.0;   // local error
	number subsetErr = 0.0;  // local contribution of marked (largest) elements

	// init iterators
	const const_iterator iterEnd = dd->template end<TElem>();
	const_iterator iter;

	// create and fill array of $\eta^2_i$ for all (local) elements
	std::vector<double> eta;
	eta.resize(numElemLocal);

	CreateListOfElemWeights<TElem>(aaError,dd->template begin<TElem>(), iterEnd, eta);
/*	size_t i=0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		const double elemErr = aaError[*iter];

		//	if no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (elemErr < 0) continue;

		eta[i++]  = elemErr;
		localErr += elemErr;
	}*/

	// sort descending using default comparison
	//std::sort (eta.begin(), eta.end(), std::greater<double>());
	UG_ASSERT(numElemLocal==eta.size(), "Huhh: number of elements does not match!");


	// error for marked subset
	subsetErr = 0.0;
	int i=0;

	// discard a fraction of elements
	UG_ASSERT( ((m_eps>=0.0) && (m_eps<=1.0)), "Huhh: m_eps invalid!");
	const int ndiscard = (m_eps!=0.0) ? (int) (numElemLocal*m_eps) : 0;
	for (i=0; i<ndiscard; ++i)
	{ subsetErr += eta[i]; }

	// define cutoff threshold
	const number theta = (m_theta*numElemLocal)/numElem;
	while (subsetErr < theta*errTotal)
	{ subsetErr += eta[i++]; }
	maxElemErr = eta[i-1];

	UG_LOG("  +++  localGoalErr^2= "<<   theta*errTotal << std::endl);
	UG_LOG("  +++  localErr^2= "<<  localErr << std::endl);
	UG_LOG("  +++  markedError^2= "<<  subsetErr << std::endl);
	UG_LOG("  +++  threshold= "<< maxElemErr<<"; i="<< i<< std::endl);


	//	mark elements with maximal contribution
	size_t numMarkedRefine = 0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		//	get element
		const double elemErr = aaError[*iter];

		// skip invalid
		if (elemErr < 0) continue;

		// skip marked
		//	if (refiner.get_mark(*iter) == RM_REFINE) continue;

		//maxElemErr = std::max(maxElemErr, elemErr);

		if (elemErr > maxElemErr)
		{
			refiner.mark(*iter, RM_REFINE);
			numMarkedRefine++;
		}
	}



#ifdef UG_PARALLEL
	if (pcl::NumProcs() > 1)
	{
		pcl::ProcessCommunicator com;
		std::size_t numMarkedRefineLocal = numMarkedRefine;
		numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
	}
#endif

	UG_LOG("  +++ EquilibrationMarkingStrategy: Marked for refinement: "<<  numMarkedRefine << " elements.\n");

}


/// marks elements above an absolute threshold (based on S. Reisters idea)
template <typename TDomain>
class AbsoluteMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	AbsoluteMarking(number eta) : m_eta(eta), m_max_level(100) {};

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

		this->m_latest_error = errTotal;

		UG_LOG("  +++ AbsoluteMarking: Error**2 = "<<errTotal <<"; marked "<< numMarkedRefine << " elements for refinement "<<std::endl);

}

}//	end of namespace

#endif




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

	virtual ~IElementMarkingStrategy(){};

	virtual void mark(elem_accessor_type& aaError,
			IRefiner& refiner,
			ConstSmartPtr<DoFDistribution> dd) = 0;
};

/// M. Breits standard refienment strategy
template <typename TDomain>
class StdMarkingStrategy : public IElementMarkingStrategy<TDomain>
{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;

	StdMarkingStrategy(number tol, number frac, int max_level)
	: m_tol(tol), m_frac(frac), m_max_level(max_level) {};

	void mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd);

protected:

	number m_tol;
	number m_frac;
	int m_max_level;
};

template <typename TDomain>
void StdMarkingStrategy<TDomain>::mark(typename base_type::elem_accessor_type& aaError,
				IRefiner& refiner,
				ConstSmartPtr<DoFDistribution> dd)
{
	MarkElementsForRefinement<typename base_type::elem_type>(aaError, refiner, dd, m_tol, m_frac, m_max_level);
}


// marks elements above a certain fraction of the maximum
template <typename TDomain>
class MaximumMarking : public IElementMarkingStrategy<TDomain>{

public:
	typedef IElementMarkingStrategy<TDomain> base_type;
	MaximumMarking(number theta=1.0) : m_theta(theta), m_eps(0.1), m_max_level(100) {};
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
		number minElemErr, maxElemErr, totalErr;
		size_t numElem;
		ComputeMinMaxTotal(aaError, dd, minElemErr, maxElemErr, totalErr, numElem);

		// init iterators
		const_iterator iter;
		const const_iterator iterEnd = dd->template end<TElem>();

		UG_LOG("  +++ Found max "<<  maxElemErr << ".\n");

		// Verfuerths strategy for skipping excess
		int discard = 0;
		const int ndiscard = (int)( numElem*m_eps );
		while (true)
		{
			// find elemnts that should be discarded
			double newmax = 0.0;
			for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
			{
				const double elemErr = aaError[*iter];

				if ( elemErr >= maxElemErr)
				{ discard++;}
				else
				{ newmax = std::max(newmax, elemErr);}
			}

			// quit, if enough elements have been found
			if ( discard <  ndiscard)
				maxElemErr = newmax;
			else
				break;
		}

		UG_LOG("  +++ Skipping " << ndiscard << " elements; new max." << maxElemErr << ".\n");

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
		}
	#endif

		UG_LOG("  +++ MaximumMarking: Marked for refinement: "<< ndiscard << "+" << numMarkedRefine << " elements.\n");

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
	number minElemErr, maxElemErr;
	number totalErr = 0.0;
	number markedErr = 0.0;  // contribution of marked (largest) elements
	size_t numElem = 0;

	ComputeMinMaxTotal(aaError, dd, minElemErr, maxElemErr, totalErr, numElem);
	//UG_LOG("  +++ Error^2= "<<  totalErr << ".\n");

	// init iterators
	const const_iterator iterEnd = dd->template end<TElem>();
	const_iterator iter;

	// TODO: Parallel
	// create array of $\eta^2_i$
	std::vector<double> eta(numElem);
	size_t i=0;
	for (iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
	{
		const double elemErr = aaError[*iter];

		//	if no error value exists: ignore (might be newly added by refinement);
		//	newly added elements are supposed to have a negative error estimator
		if (elemErr < 0) continue;

		eta[i++] = elemErr;
		// totalErr += elemErr;
	}

	UG_ASSERT(i==numElem, "Huhh: number of elements does not match!");
	UG_ASSERT( ((m_eps>=0.0) && (m_eps<=1.0)), "Huhh: m_eps invalid!");

	// sort descending using default comparison
	std::sort (eta.begin(), eta.end(), std::greater<double>());

	// discard at least some fraction of elements
	const int ndiscard = (m_eps!=0.0) ? (int) (numElem*m_eps) : 0;

	// find cutoff threshold
	UG_LOG("  +++  Error^2= "<<   m_theta*totalErr << std::endl);
	markedErr = 0.0;
	i=0;
	while ( markedErr < m_theta*totalErr)
	{
		markedErr += eta[i++];
	}
	maxElemErr = eta[i-1];

	UG_LOG("  +++  Error^2= "<<  totalErr << "; MarkedError^2= "<<  markedErr << "; threshold= "<< maxElemErr<<"; i="<< i<<"\n");


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
			markedErr += elemErr;
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

		size_t numMarkedRefine = 0;
		for (const_iterator iter = dd->template begin<TElem>(); iter != iterEnd; ++iter)
		{
		//	get element error
			TElem* elem = *iter;
			const number elemError = aaError[elem];

		// skip newly added
			if (elemError < 0) continue;

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
			numMarkedRefine = com.allreduce(numMarkedRefineLocal, PCL_RO_SUM);
		}
	#endif

		UG_LOG("  +++ AbsoluteMarking: Marked "<< numMarkedRefine << " elements for refinement.\n");

}

}//	end of namespace

#endif


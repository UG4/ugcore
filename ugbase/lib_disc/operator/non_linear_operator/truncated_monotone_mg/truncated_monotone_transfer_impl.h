/*
 * truncated_monotone_transfer_impl.h
 *
 *  Created on: 20.01.2014
 *      Author: raphaelprohl
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__TRUNCATED_MONOTONE_MG__TRUNCATED_MONOTONE_TRANSFER_IMPL_H_
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__TRUNCATED_MONOTONE_MG__TRUNCATED_MONOTONE_TRANSFER_IMPL_H_

#include "truncated_monotone_transfer.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
void
TruncatedMonotoneTransfer<TDomain, TAlgebra>::
set_levels(GridLevel coarseLevel, GridLevel fineLevel)
{
	m_bInit = false;
	m_fineLevel = fineLevel;
	m_coarseLevel = coarseLevel;

	if(m_fineLevel.level () - m_coarseLevel.level () != 1)
		UG_THROW("TruncatedMonotoneTransfer: Can only project between successive level.");
}

template <typename TDomain, typename TAlgebra>
void
TruncatedMonotoneTransfer<TDomain, TAlgebra>::
init()
{

	//	Done:
	m_bInit = true;
}

template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
TruncatedMonotoneTransfer<TDomain, TAlgebra>::
prolongation(const GridLevel& fineGL, const GridLevel& coarseGL,
			 ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
	SmartPtr<matrix_type> stdP = base_type::prolongation(fineGL, coarseGL, spApproxSpace);

	return stdP;
}

template <typename TDomain, typename TAlgebra>
SmartPtr<typename TAlgebra::matrix_type>
TruncatedMonotoneTransfer<TDomain, TAlgebra>::
restriction(const GridLevel& coarseGL, const GridLevel& fineGL,
		    ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace)
{
	SmartPtr<matrix_type> stdR = base_type::restriction(coarseGL, fineGL, spApproxSpace);

	//TODO: modify stdR for activeDofs! IConstraint::adjust_restriction?
	//		or only for the special case IObstacleConstraint?

	return stdR;
}

template <typename TDomain, typename TAlgebra>
SmartPtr<ITransferOperator<TDomain, TAlgebra> >
TruncatedMonotoneTransfer<TDomain, TAlgebra>::
clone()
{
	SmartPtr<TruncatedMonotoneTransfer> op (new TruncatedMonotoneTransfer());
	for(size_t i = 0; i < m_vConstraint.size (); ++i)
		op->add_constraint (m_vConstraint[i]);

	//TODO: should anything else be done for cloning?!
	return op;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__TRUNCATED_MONOTONE_MG__TRUNCATED_MONOTONE_TRANSFER_IMPL_H_ */

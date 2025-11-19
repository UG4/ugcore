/*
 * Copyright (c) 2014-2015:  G-CSC, Goethe University Frankfurt
 * Author: Raphael Prohl
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
	//		(only for the "finest grid to (finest-1) grid-transfer; then: canonical)

	int topLev = spApproxSpace->num_levels() - 1;
	if(topLev == fineGL.level())
	{
		UG_LOG("on topLevel: " <<topLev<<"\n");

		ConstSmartPtr<DoFDistribution> spCoarseDD = spApproxSpace->dof_distribution(coarseGL);
		ConstSmartPtr<DoFDistribution> spFineDD = spApproxSpace->dof_distribution(fineGL);

		#ifdef UG_PARALLEL
		stdR->set_storage_type(PST_CONSISTENT);
		#endif

		// 	adjust using constraints
		for(size_t i = 0; i < m_vConstraint.size(); ++i)
			if(m_vConstraint[i]->type() == CT_CONSTRAINTS)
			{
				UG_LOG("is CT_CONSTRAINT")
				m_vConstraint[i]->adjust_restriction(*stdR, spCoarseDD, spFineDD, CT_CONSTRAINTS);
			}
	}

	//TODO: restrict the obstacle-value in a monotone way! (in every level-transfer)
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

#endif
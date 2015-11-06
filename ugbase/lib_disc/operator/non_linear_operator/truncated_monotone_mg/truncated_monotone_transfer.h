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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__TRUNCATED_MONOTONE_MG__TRUNCATED_MONOTONE_TRANSFER_H_
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__TRUNCATED_MONOTONE_MG__TRUNCATED_MONOTONE_TRANSFER_H_

#include "lib_disc/operator/linear_operator/std_transfer.h"

namespace ug{

template <typename TDomain, typename TAlgebra>
class TruncatedMonotoneTransfer:
	public StdTransfer<TDomain, TAlgebra>
{
	public:
	/// This type
		typedef TruncatedMonotoneTransfer<TDomain, TAlgebra> this_type;

	/// base type
		typedef StdTransfer<TDomain, TAlgebra> base_type;

	///	Type of Domain
		typedef TDomain domain_type;

	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Matrix
		typedef typename TAlgebra::matrix_type matrix_type;

	public:
	///	Constructor
		TruncatedMonotoneTransfer(): ITransferOperator<TDomain, TAlgebra>(),
				m_bInit(false)
		{};

	public:
		//////////////////////////////
		//	INTERFACE METHODS
		//////////////////////////////

	///	initialize the operator
		void init();

	/// Set Levels for Prolongation coarse -> fine
		void set_levels(GridLevel coarseLevel, GridLevel fineLevel);

	///	returns prolongation as a matrix
		SmartPtr<matrix_type> prolongation(const GridLevel& fineGL,
				const GridLevel& coarseGL,
				ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace);

	///	returns restriction as a matrix
		SmartPtr<matrix_type> restriction(const GridLevel& coarseGL,
				const GridLevel& fineGL,
				ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace);
	///	Clone
		SmartPtr<ITransferOperator<TDomain, TAlgebra> > clone();

	private:
	///	init flag
		bool m_bInit;

	///	list of post processes
		using base_type::m_vConstraint;

	///	coarse and fine Gridlevel
		GridLevel m_coarseLevel;
		GridLevel m_fineLevel;
};

} // end namespace ug

// include implementation
#include "truncated_monotone_transfer_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__TRUNCATED_MONOTONE_MG__TRUNCATED_MONOTONE_TRANSFER_H_ */

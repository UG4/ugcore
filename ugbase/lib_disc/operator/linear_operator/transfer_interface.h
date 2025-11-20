/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_INTERFACE__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_INTERFACE__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_grid/tools/grid_level.h"
#include "lib_disc/function_spaces/grid_function.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Transfer Operator
///////////////////////////////////////////////////////////////////////////////

/// interface for transfer routines
template <typename TDomain, typename TAlgebra>
class ITransferOperator
{
	public:
	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	///	Domain type
		using domain_type = TDomain;

	public:
	///	constructor
		ITransferOperator(){clear_constraints();}

	///	clears dirichlet post processes
		virtual void clear_constraints(){m_vConstraint.clear();};

	///	adds a dirichlet post process (not added if already registered)
		virtual void add_constraint(SmartPtr<IConstraint<TAlgebra> > pp){
			//	add only once
			if(std::find(m_vConstraint.begin(), m_vConstraint.end(), pp) !=
					m_vConstraint.end()) return;
			m_vConstraint.push_back(pp);
		};

	///	removes a post process
		virtual void remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp){
			m_vConstraint.erase(m_vConstraint.begin(),
			    std::remove(m_vConstraint.begin(), m_vConstraint.end(), pp));
		}

	public:
	///	initialize the operator
		virtual void init() = 0;

	/// Set Levels for Prolongation coarse -> fine
		virtual void set_levels(GridLevel coarseLevel, GridLevel fineLevel) = 0;

	/// Prolongates vector, i.e. moves data from coarse to fine level
		virtual void prolongate(vector_type& uFine, const vector_type& uCoarse) = 0;

	/// Restricts vector, i.e. moves data from fine to coarse level
		virtual void do_restrict(vector_type& uCoarse, const vector_type& uFine) = 0;

	///	returns prolongation as a matrix
		virtual SmartPtr<matrix_type>
		prolongation(const GridLevel& fineGL, const GridLevel& coarseGL,
		             ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace){
			UG_THROW("ITransferOperator: Matrix-prolongation not implemented.")
		}

	///	returns restriction as a matrix
		virtual SmartPtr<matrix_type>
		restriction(const GridLevel& coarseGL, const GridLevel& fineGL,
		            ConstSmartPtr<ApproximationSpace<TDomain> > spApproxSpace){
			UG_THROW("ITransferOperator: Matrix-restriction not implemented.")
		}

	///	Clone
		virtual SmartPtr<ITransferOperator > clone() = 0;

	///	virtual destructor
		virtual ~ITransferOperator() = default;

	protected:
	///	list of post processes
		std::vector<SmartPtr<IConstraint<TAlgebra> > > m_vConstraint;

};

///////////////////////////////////////////////////////////////////////////////
// Transfer Post Process
///////////////////////////////////////////////////////////////////////////////

/// interface for transfer routines
template <typename TDomain, typename TAlgebra>
class ITransferPostProcess
{
	public:
	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Domain type
		using domain_type = TDomain;

	///	GridFunction type
		using GF = GridFunction<TDomain, TAlgebra>;

	public:
	/// apply post process
		virtual void post_process(SmartPtr<GF> spGF) = 0;

	///	virtual destructor
		virtual ~ITransferPostProcess() = default;
};

} // end namespace ug

#endif
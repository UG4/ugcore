/*
 * transfer_interface.h
 *
 *  Created on: 20.12.2011
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_INTERFACE__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_INTERFACE__

#include "lib_algebra/operator/interface/operator.h"
#include "lib_disc/dof_manager/grid_level.h"

namespace ug{

//predeclaration
template <typename TAlgebra>
class IConstraint;


///////////////////////////////////////////////////////////////////////////////
// Transfer Operator
///////////////////////////////////////////////////////////////////////////////

template <typename TAlgebra>
class ITransferOperator
{
	public:
	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	public:
	/// Set Levels for Prolongation coarse -> fine
		virtual void set_levels(GridLevel coarseLevel, GridLevel fineLevel) = 0;

	///	clears dirichlet post processes
		virtual void clear_constraints() = 0;

	///	adds a dirichlet post process (not added if already registered)
		virtual void add_constraint(SmartPtr<IConstraint<TAlgebra> > pp) = 0;

	///	removes a post process
		virtual void remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp) = 0;

	public:
	///	initialize the operator
		virtual void init() = 0;

	/// Apply Operator, interpolate function, prolongate
		virtual void apply(vector_type& uFine, const vector_type& uCoarse) = 0;

	/// Apply Transposed Operator u = L^T*f, restrict
		virtual void apply_transposed(vector_type& uCoarse, const vector_type& uFine) = 0;

	///	Clone
		virtual SmartPtr<ITransferOperator<TAlgebra> > clone() = 0;

	///	virtual destructor
		~ITransferOperator() {}
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__TRANSFER_INTERFACE__ */

/*
 * projection_operator.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "lib_grid/lg_base.h"
#include "transfer_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

/**
 * The Projection operator transfers is used to transfer vectors between two
 * grid levels. It implements a purely algebraic interface, just mapping
 * between two algebraic vectors, but given the approximation space this is indeed
 * a mapping between two grid functions.
 *
 * \tparam	TDomain		the domain
 * \tparam	TAlgebra	the algebra
 */
template <typename TDomain, typename TAlgebra>
class InjectionTransfer :
	virtual public ITransferOperator<TAlgebra>
{
	public:
	///	Type of algebra
		typedef TAlgebra algebra_type;

	///	Type of Vector
		typedef typename TAlgebra::vector_type vector_type;

	///	Type of Vector
		typedef typename TAlgebra::matrix_type matrix_type;

	///	Type of Domain
		typedef TDomain domain_type;

	public:
	///	Constructor
		InjectionTransfer() : m_bInit(false) {}

	///	Constructor
		InjectionTransfer(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
			m_spApproxSpace(approxSpace), m_bInit(false)
		{}

	///	Set Approximation Space
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	///	virtual Destructor
		virtual ~InjectionTransfer(){};
	public:
	///	Set approximation level
		void set_levels(GridLevel coarseLevel, GridLevel fineLevel);

	///	clears dirichlet post processes
		void clear_constraints() {m_vConstraint.clear();}

	///	adds a dirichlet post process (not added if already registered)
		void add_constraint(SmartPtr<IConstraint<TAlgebra> > pp);

	///	removes a post process
		void remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp);

	public:
	///	Init operator
		virtual void init();

	/// Project uFine to uCoarse; uCoarse = P(uFine);
		virtual void prolongate(vector_type& uFine, const vector_type& uCoarse);

	/// Apply Transposed Operator u = L^T*f
		virtual void restrict(vector_type& uCoarse, const vector_type& uFine);

	///	clones the operator
		virtual SmartPtr<ITransferOperator<TAlgebra> > clone();

	protected:
	/// matrix used for projection
		matrix_type m_matrix;

	///	the underlying approximation space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	list of post processes
		std::vector<SmartPtr<IConstraint<TAlgebra> > > m_vConstraint;

	///	fine level of approximation space
		GridLevel m_fineLevel;

	///	coarse level of approximation space
		GridLevel m_coarseLevel;

	///	init flag
		bool m_bInit;
};

}

#include "projection_operator_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__ */

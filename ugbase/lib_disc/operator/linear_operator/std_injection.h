/*
 * std_injection.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_INJECTION__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_INJECTION__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "transfer_interface.h"

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
class StdInjection :
	virtual public ITransferOperator<TDomain, TAlgebra>
{
	public:
	///	Type of base class
		typedef ITransferOperator<TDomain, TAlgebra> base_type;

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
		StdInjection() : m_bInit(false) {}

	///	Constructor
		StdInjection(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
			m_spApproxSpace(approxSpace), m_bInit(false)
		{}

	///	Set Approximation Space
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	///	virtual Destructor
		virtual ~StdInjection(){};
	public:
	///	Set approximation level
		void set_levels(GridLevel coarseLevel, GridLevel fineLevel);

	protected:
		template <typename TElem>
		void set_identity_on_pure_surface(matrix_type& mat,
		                                  const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

		void set_identity_on_pure_surface(matrix_type& mat,
		                                  const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

	public:
	///	Init operator
		virtual void init();

	/// Project uFine to uCoarse; uCoarse = P(uFine);
		virtual void prolongate(vector_type& uFine, const vector_type& uCoarse);

	/// Apply Transposed Operator u = L^T*f
		virtual void do_restrict(vector_type& uCoarse, const vector_type& uFine);

	///	clones the operator
		virtual SmartPtr<ITransferOperator<TDomain, TAlgebra> > clone();

	protected:
	/// matrix used for projection
		matrix_type m_matrix;

	///	the underlying approximation space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	fine level of approximation space
		GridLevel m_fineLevel;

	///	coarse level of approximation space
		GridLevel m_coarseLevel;

	///	init flag
		bool m_bInit;
};

}

#include "std_injection_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__STD_INJECTION__ */

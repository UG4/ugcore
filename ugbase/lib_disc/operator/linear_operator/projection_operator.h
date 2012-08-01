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
 * This functions assembles the interpolation matrix between to
 * grid levels using only the Vertex degrees of freedom.
 *
 * \param[out]	mat 			Assembled interpolation matrix that interpolates u -> v
 * \param[in] 	approxSpace		Approximation Space
 * \param[in]	coarseLevel		Coarse Level index
 * \param[in]	fineLevel		Fine Level index
 */
template <typename TDD, typename TAlgebra>
void AssembleVertexProjection(typename TAlgebra::matrix_type& mat,
                              const TDD& coarseDD, const TDD& fineDD);

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
class P1Projection :
	virtual public IProjectionOperator<	typename TAlgebra::vector_type,
										typename TAlgebra::vector_type>
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
		P1Projection() : m_bInit(false) {}

	///	Constructor
		P1Projection(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
			m_spApproxSpace(approxSpace), m_bInit(false)
		{}

	///	Set Approximation Space
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	///	Set approximation level
		void set_levels(GridLevel coarseLevel, GridLevel fineLevel);

	public:
	///	Init operator
		virtual void init();

		virtual void init(const vector_type& u){init();}

	/// Project uFine to uCoarse; uCoarse = P(uFine);
		virtual void apply(vector_type& uCoarseOut, const vector_type& uFineIn);

	/// Apply Transposed Operator u = L^T*f
		virtual void apply_transposed(vector_type& uFineOut, const vector_type& uCoarseIn);

	/// Apply sub not implemented
		virtual void apply_sub(vector_type& u, const vector_type& v);

	///	clones the operator
		virtual SmartPtr<IProjectionOperator<vector_type> > clone();

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

#include "projection_operator_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROJECTION_OPERATOR__ */

/*
 * prolongation_operator.h
 *
 *  Created on: 04.12.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__

// extern headers
#include <iostream>

// other ug4 modules
#include "common/common.h"
#include "transfer_interface.h"
#include "lib_disc/spatial_disc/constraints/constraint_interface.h"
#include "lib_algebra/operator/debug_writer.h"

#ifdef UG_PARALLEL
#include "lib_disc/parallelization/parallelization_util.h"
#endif

namespace ug{

///	Prologation Operator for P1 Approximation Spaces
/**	By default a special optimization is performed for p1-lagrange-elements.
 * This optimization is only valid if all elements have been refined with
 * standard refinement rules. If closure elements are generated, this optimization
 * has to be deactivated (use StdTransfer::enable_p1_lagrange_optimization(false)).
 */
template <typename TDomain, typename TAlgebra>
class StdTransfer :
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
	/// Default constructor
		StdTransfer() : m_bInit(false), m_p1LagrangeOptimizationEnabled(true),
						m_dampRes(1.0), m_spDebugWriter(NULL)
		{clear_constraints();};

	///	Constructor setting approximation space
		StdTransfer(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
			m_spApproxSpace(approxSpace), m_bInit(false),
			m_p1LagrangeOptimizationEnabled(true), m_dampRes(1.0), m_spDebugWriter(NULL)
		{clear_constraints();};

	///	Set approximation space
		void set_approximation_space(SmartPtr<ApproximationSpace<TDomain> > approxSpace);

	///	set interpolation damping
		void set_restriction_damping(number damp) {m_dampRes = damp;}

	///	set debug writer
		void set_debug(SmartPtr<IDebugWriter<TAlgebra> > spDebugWriter) {
			m_spDebugWriter = spDebugWriter;
		}

	/// virtual destructor
		virtual ~StdTransfer(){};

	public:
	///	Set levels
		virtual void set_levels(GridLevel coarseLevel, GridLevel fineLevel);

	///	clears dirichlet post processes
		void clear_constraints() {m_vConstraint.clear();}

	///	adds a dirichlet post process (not added if already registered)
		void add_constraint(SmartPtr<IConstraint<TAlgebra> > pp);

	///	removes a post process
		void remove_constraint(SmartPtr<IConstraint<TAlgebra> > pp);

	///	enables/disables an assembling optimization for p1-lagrange elements
	/**	The optimization is enabled by default. It can however only be used,
	 * if all elements are refined with their standard refinement rule. If one
	 * uses anisotropic refinement or refinement with closure, the optimization
	 * should be disabled.
	 * \todo	The normal assembling strategy should be optimized in such a way
	 * 			that the p1-lagrange optimization is no longer required. This
	 * 			however involves something like a ref-type-hash for each element,
	 * 			which returns a unique number based on the types and order of children.*/
		void enable_p1_lagrange_optimization(bool enable)	{m_p1LagrangeOptimizationEnabled = enable;}
		bool p1_lagrange_optimization_enabled() const		{return m_p1LagrangeOptimizationEnabled;}

	public:
	///	initialize the operator
		virtual void init();

	/// apply Operator, interpolate function
		virtual void prolongate(vector_type& uFineOut, const vector_type& uCoarse);

	/// apply transposed Operator, restrict function
		virtual void do_restrict(vector_type& uCoarse, const vector_type& uFine);

	///	returns new instance with same setting
		virtual SmartPtr<ITransferOperator<TAlgebra> > clone();

	protected:
	///	debug writing of matrix
		void write_debug(const matrix_type& mat, const char* filename);

		void assemble_restriction_elemwise(matrix_type& mat,
		                                    const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
		                                    ConstSmartPtr<TDomain> spDomain);

		void assemble_restriction_p1(matrix_type& mat,
		                              const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

		template <typename TElem>
		void set_identity_on_pure_surface(matrix_type& mat,
		                                  const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

		void set_identity_on_pure_surface(matrix_type& mat,
		                                  const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

	protected:
	///	matrix to store prolongation
		matrix_type m_Restriction;

	///	list of post processes
		std::vector<SmartPtr<IConstraint<TAlgebra> > > m_vConstraint;

	///	approximation space
		SmartPtr<ApproximationSpace<TDomain> > m_spApproxSpace;

	///	fine grid level
		GridLevel m_fineLevel;

	///	coarse grid level
		GridLevel m_coarseLevel;

	///	initialization flag
		bool m_bInit;

	///	flag for p1-lagrange-optimization
		bool m_p1LagrangeOptimizationEnabled;

	///	damping parameter
		number m_dampRes;

	///	debug writer
		SmartPtr<IDebugWriter<TAlgebra> > m_spDebugWriter;
};

} // end namespace ug

#include "prolongation_operator_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__ */

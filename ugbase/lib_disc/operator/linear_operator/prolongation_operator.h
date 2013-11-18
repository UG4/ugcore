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
		StdTransfer() : m_bInit(false), m_dampRes(1.0), m_spDebugWriter(NULL) {clear_constraints();};

	///	Constructor setting approximation space
		StdTransfer(SmartPtr<ApproximationSpace<TDomain> > approxSpace) :
			m_spApproxSpace(approxSpace), m_bInit(false), m_dampRes(1.0), m_spDebugWriter(NULL)
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

		void assemble_prolongation_elemwise(typename TAlgebra::matrix_type& mat,
		                                    const DoFDistribution& coarseDD, const DoFDistribution& fineDD,
		                                    ConstSmartPtr<TDomain> spDomain);

		void assemble_prolongation_p1(typename TAlgebra::matrix_type& mat,
		                              const DoFDistribution& coarseDD, const DoFDistribution& fineDD);

	protected:
	///	matrix to store prolongation
		matrix_type m_matrix;

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

	///	damping parameter
		number m_dampRes;

	///	debug writer
		SmartPtr<IDebugWriter<TAlgebra> > m_spDebugWriter;
};

} // end namespace ug

#include "prolongation_operator_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__PROLONGATION_OPERATOR__ */

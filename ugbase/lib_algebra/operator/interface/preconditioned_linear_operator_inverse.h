/*
 * preconditioned_linear_operator_inverse.h
 *
 *  Created on: 23.04.2013
 *      Author: andreasvogel
 */

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONED_LINEAR_OPERATOR_INVERSE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONED_LINEAR_OPERATOR_INVERSE__

#include "linear_operator_inverse.h"
#include "linear_iterator.h"
#include "lib_algebra/operator/debug_writer.h"
#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/convergence_check.h"
#include "common/log.h"
#include "linear_solver_profiling.h"

namespace ug{
///////////////////////////////////////////////////////////////////////////////
// Inverse of a Linear Operator using a ILinearIterator as preconditioner
///////////////////////////////////////////////////////////////////////////////

/// describes an inverse linear mapping X->X
/**
 * This a useful derived class from ILinearOperatorInverse, that uses a
 * ILinearIterator in order to precondition the solution process. This is
 * used e.g. in LinearSolver, CG and BiCGStab.
 *
 * \tparam	X		domain and range space
 */
template <typename X>
class IPreconditionedLinearOperatorInverse
	: public ILinearOperatorInverse<X>,
	  public VectorDebugWritingObject<X>
{
	public:
	///	Domain space
		typedef X domain_function_type;

	///	Range space
		typedef X codomain_function_type;

	///	Base class
		typedef ILinearOperatorInverse<X,X> base_type;

	protected:
		using base_type::linear_operator;

	public:
		using VectorDebugWritingObject<X>::write_debug;
		using base_type::name;
		using base_type::apply_return_defect;

	public:
	///	Empty constructor
		IPreconditionedLinearOperatorInverse()
			: m_bRecompute(false), m_spPrecond(NULL)
		{}

	///	constructor setting the preconditioner
		IPreconditionedLinearOperatorInverse(SmartPtr<ILinearIterator<X,X> > spPrecond)
			: m_bRecompute(false), m_spPrecond(spPrecond)
		{}

	///	constructor setting the preconditioner
		IPreconditionedLinearOperatorInverse(SmartPtr<ILinearIterator<X,X> > spPrecond,
		                                     SmartPtr<IConvergenceCheck<X> > spConvCheck)
			: 	base_type(spConvCheck),
				m_bRecompute(false), m_spPrecond(spPrecond)
		{}

	///	sets the preconditioner
		void set_preconditioner(SmartPtr<ILinearIterator<X, X> > spPrecond)
		{
			m_spPrecond = spPrecond;
		}

	///	returns the preconditioner
		SmartPtr<ILinearIterator<X, X> > preconditioner(){return m_spPrecond;}

	///	initializes the solver for an operator
		virtual bool init(SmartPtr<ILinearOperator<X,X> > J, const X& u)
		{
			if(!base_type::init(J, u)) return false;

			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(J, u))
					UG_THROW(name() << "::init: Cannot init Preconditioner "
													"Operator for Operator J.");
			LS_PROFILE_END();

			return true;
		}

	///	initializes the solver for an operator
		virtual bool init(SmartPtr<ILinearOperator<X,X> > L)
		{
			if(!base_type::init(L)) return false;

			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(L))
					UG_THROW(name() <<"::prepare: Cannot init Preconditioner "
														"Operator for Operator L.");
			LS_PROFILE_END();

			return true;
		}

		virtual bool apply(X& x, const X& b)
		{
		//	copy defect
			SmartPtr<X> spB = b.clone(); X& bTmp = *spB;
//			X bTmp; bTmp.resize(b.size()); bTmp = b;

		//	solve on copy of defect
			bool bRes = apply_return_defect(x, bTmp);

		//	write updated defect
			write_debug(bTmp, "LS_UpdatedDefectEnd.vec");

		//	compute defect again, for debug purpose
			if(m_bRecompute)
			{
			//	recompute defect
				bTmp = b; linear_operator()->apply_sub(bTmp, x);
				number norm = bTmp.norm();

			//	print norm of recomputed defect
				UG_LOG("%%%% DEBUG "<<name()<<": (Re)computed defect has norm: "
				       <<norm<<"\n");

			//	write true end defect
				write_debug(bTmp, "LS_TrueDefectEnd.vec");
			}

		//	return
			return bRes;
		}

	///	for debug: computes norm again after whole calculation of apply
		void set_compute_fresh_defect_when_finished(bool bRecompute)
		{
			m_bRecompute = bRecompute;
		}

	protected:
	///	flag if fresh defect should be computed when finish for debug purpose
		bool m_bRecompute;

	///	Iterator used in the iterative scheme to compute the correction and update the defect
		SmartPtr<ILinearIterator<X,X> > m_spPrecond;

};

}
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONED_LINEAR_OPERATOR_INVERSE__ */

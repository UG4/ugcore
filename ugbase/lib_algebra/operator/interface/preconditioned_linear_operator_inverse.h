/*
 * Copyright (c) 2013-2015:  G-CSC, Goethe University Frankfurt
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

#ifndef __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONED_LINEAR_OPERATOR_INVERSE__
#define __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONED_LINEAR_OPERATOR_INVERSE__

#include "linear_operator_inverse.h"
#include "linear_iterator.h"
#include "lib_algebra/operator/debug_writer.h"
#include "common/util/smart_pointer.h"
#include "lib_algebra/operator/convergence_check.h"
#include "common/log.h"
#include "linear_solver_profiling.h"

#undef DEBUG_FOR_AMG

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
		using base_type::m_spConvCheck;

	public:
		using VectorDebugWritingObject<X>::write_debug;
		using base_type::name;
		using base_type::apply_return_defect;

	public:
	///	Empty constructor
		IPreconditionedLinearOperatorInverse()
			: m_bRecompute(false), m_spPrecond(NULL)
#ifdef DEBUG_FOR_AMG
, m_amgDebug(0)
#endif
		{}

	///	constructor setting the preconditioner
		IPreconditionedLinearOperatorInverse(SmartPtr<ILinearIterator<X,X> > spPrecond)
			: m_bRecompute(false), m_spPrecond(spPrecond)
#ifdef DEBUG_FOR_AMG
, m_amgDebug(0)
#endif
		{}

	///	constructor setting the preconditioner
		IPreconditionedLinearOperatorInverse(SmartPtr<ILinearIterator<X,X> > spPrecond,
		                                     SmartPtr<IConvergenceCheck<X> > spConvCheck)
			: 	base_type(spConvCheck),
				m_bRecompute(false), m_spPrecond(spPrecond)
#ifdef DEBUG_FOR_AMG
, m_amgDebug(0)
#endif
		{}

	///	sets the preconditioner
		void set_preconditioner(SmartPtr<ILinearIterator<X, X> > spPrecond)
		{
			m_spPrecond = spPrecond;
		}

	///	returns the preconditioner
	/// \{
		SmartPtr<ILinearIterator<X, X> > preconditioner(){return m_spPrecond;}
		ConstSmartPtr<ILinearIterator<X, X> > preconditioner() const {return m_spPrecond;}
	/// \}

	///	initializes the solver for an operator
		virtual bool init(SmartPtr<ILinearOperator<X,X> > J, const X& u)
		{
			if(!base_type::init(J, u)) return false;

			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(J, u))
					UG_THROW(name() << "::init: Cannot init Preconditioner "
													"Operator for Operator J.");
			LS_PROFILE_END(LS_InitPrecond);

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
			LS_PROFILE_END(LS_InitPrecond);

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

#ifdef DEBUG_FOR_AMG
			if (m_amgDebug>0)
			{
			// convergence post-check
			X myError(x.size());
			myError.set_random(-1.0, 1.0);

			bTmp.set(0.0);

			this->write_debug(myError, "AMGDebugPre");
			apply_return_defect(myError, bTmp);
			this->write_debug(myError, "AMGDebugPost");
			}
#endif


		//	return
			return bRes;
		}

	///	returns config information of convergence check and preconditioner
		std::string config_string_preconditioner_convergence_check() const
		{
			std::stringstream ss;
			ss << " Convergence Check: ";
			if(m_spConvCheck.valid()) ss << ConfigShift(m_spConvCheck->config_string()) << "\n";
			else ss << "  NOT SET!\n";
			ss << " Preconditioner: ";
			if(m_spPrecond.valid()) ss << ConfigShift(m_spPrecond->config_string()) << "\n";
			else ss << "  NOT SET!\n";
			return ss.str();
		}

	///	returns information about configuration parameters
		virtual std::string config_string() const
		{
			std::stringstream ss;
			ss << name() << "\n" << config_string_preconditioner_convergence_check();
			return ss.str();
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

#ifdef DEBUG_FOR_AMG
	public:
		void set_debug_amg(int b) {m_amgDebug = b;}
		int m_amgDebug;
#endif
};

}
#endif /* __H__LIB_ALGEBRA__OPERATOR__INTERFACE__PRECONDITIONED_LINEAR_OPERATOR_INVERSE__ */

/*
 * linear_solver.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/operator_interface.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

#define PROFILE_LS
#ifdef PROFILE_LS
	#define LS_PROFILE_FUNC()		PROFILE_FUNC()
	#define LS_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define LS_PROFILE_END()		PROFILE_END()
#else
	#define LS_PROFILE_FUNC()
	#define LS_PROFILE_BEGIN(name)
	#define LS_PROFILE_END()
#endif

namespace ug{

/// linear solver using abstract preconditioner interface
/**
 * This class is a linear iterating scheme, that uses any implementation
 * of the ILinearIterator interface to precondition the iteration.
 *
 * \tparam 		TAlgebra		algebra type
 */
template <typename TAlgebra>
class LinearSolver
	: public ILinearOperatorInverse<  typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	  	  typename TAlgebra::vector_type>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Base type
		typedef ILinearOperatorInverse<vector_type,vector_type> base_type;

	protected:
		using base_type::convergence_check;

	public:
	///	Default constructor
		LinearSolver() :
			m_A(NULL), m_spPrecond(NULL),
			m_bRecomputeDefectWhenFinished(false),
			m_pDebugWriter(NULL)
		{}

	///	returns the name of the solver
		virtual const char* name() const {return "Iterative Linear Solver";}

	///	for debug: computes norm again after whole calculation of apply
		void set_compute_fresh_defect_when_finished(bool bRecomputeDefectWhenFinished)
		{
			m_bRecomputeDefectWhenFinished = bRecomputeDefectWhenFinished;
		}

	///	sets the preconditioner
		void set_preconditioner(SmartPtr<ILinearIterator<vector_type, vector_type> > precond)
		{
			m_spPrecond = precond;
		}

	///	initializes the solver for an operator
		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
		//	remember operator
			m_A = &J;

		//	init the preconditioner
			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(J, u))
					UG_THROW_FATAL("LinearSolver::prepare: Cannot init Iterator "
									"Operator for Operator J.");

			LS_PROFILE_END();

		//	we're done
			return true;
		}

	///	initializes the solver for an operator
		virtual bool init(ILinearOperator<vector_type, vector_type>& L)
		{
		//	remember operator
			m_A = &L;

		// 	init Preconditioner for operator A
			LS_PROFILE_BEGIN(LS_InitPrecond);
			if(m_spPrecond.valid())
				if(!m_spPrecond->init(L))
					UG_THROW_FATAL("LinearSolver::prepare: Cannot init Iterator "
									"Operator for Operator L.");
			LS_PROFILE_END();

		//	we're done
			return true;
		}

	///	solves the system and returns the last defect
		virtual bool apply_return_defect(vector_type& x, vector_type& b)
		{
			LS_PROFILE_BEGIN(LS_ApplyReturnDefect);

			if(m_A == NULL)
			{
				UG_LOG("ERROR in 'LinearSolver::apply': "
						"Operator that should be inverted has not been set.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				UG_THROW_FATAL("LinearSolver::apply: Inadequate storage format of Vectors.");
			#endif

		// 	rename b as d (for convenience)
			vector_type& d = b;

		// 	build defect:  d := b - J(u)*x
			LS_PROFILE_BEGIN(LS_BuildDefect);
			m_A->apply_sub(d, x);
			LS_PROFILE_END(); //LS_BuildDefect

		// 	create correction
		// 	todo: 	it would be sufficient to only copy the pattern (and parallel constructor)
		//			without initializing the values
			LS_PROFILE_BEGIN(LS_CreateCorrection);
			vector_type c; c.create(x.size()); c = x;
			LS_PROFILE_END();

			LS_PROFILE_BEGIN(LS_ComputeStartDefect);
			prepare_conv_check();
			convergence_check()->start(d);
			LS_PROFILE_END();

		// 	Iteration loop
			while(!convergence_check()->iteration_ended())
			{
			// 	Compute a correction c := B*d using one iterative step
			// 	Internally the defect is updated d := d - A*c = d - A*(x+c)
				if(m_spPrecond.valid()) {
					LS_PROFILE_BEGIN(LS_ApplyPrecond);
					if(!m_spPrecond->apply_update_defect(c, d))
					{
						UG_LOG("ERROR in 'LinearSolver::apply': Iterator "
								"Operator applied incorrectly. Aborting.\n");
						return false;
					}
					LS_PROFILE_END(); //LS_ApplyPrecond
				}
				write_debug_vector(d, "Defect.vec");
				write_debug_vector(c, "Correction.vec");

			// 	add correction to solution: x += c
				LS_PROFILE_BEGIN(LS_AddCorrection);
				x += c;
				LS_PROFILE_END(); //LS_AddCorrection

			// 	compute new defect (in parallel) d := d - A*c
				LS_PROFILE_BEGIN(LS_ComputeNewDefect);
				convergence_check()->update(d);
				LS_PROFILE_END(); //LS_ComputeNewDefect
			}

		//	write some information when ending the iteration
			if(!convergence_check()->post())
			{
				UG_LOG("ERROR in 'LinearSolver::apply': post-convergence-check "
						"signaled failure. Aborting.\n");
				return false;
			}

		//	end profiling of whole function
			LS_PROFILE_END(); //LS_ApplyReturnDefect

		//	we're done
			return true;
		}

		virtual bool apply(vector_type& x, const vector_type& b)
		{
		//	copy defect
			vector_type b2; b2.resize(b.size());
			b2 = b;

		//	solve on copy of defect
			bool bRes = apply_return_defect(x, b2);

			write_debug_vector(b2, "LS_UpdatedDefectEnd.vec");

		//	compute defect again, for debug purpose
			if(m_bRecomputeDefectWhenFinished)
			{
				b2 = b;
				m_A->apply_sub(b2, x);

				number norm = b2.two_norm();
				UG_LOG("%%%% DEBUG: (Re)computed defect has norm: "<<norm<<"\n");

				write_debug_vector(b2, "LS_TrueDefectEnd.vec");
			}

		//	return
			return bRes;
		}

		// destructor
		virtual ~LinearSolver() {};

	///	set debug writer
		void set_debug(IDebugWriter<algebra_type>* debugWriter)
		{
			m_pDebugWriter = debugWriter;
		}

	///	returns the debug writer
		IDebugWriter<algebra_type>* debug_writer() {return m_pDebugWriter;}

	protected:
	///	prepares the convergence check output
		void prepare_conv_check()
		{
			convergence_check()->set_name(name());
			convergence_check()->set_symbol('%');
			if(m_spPrecond.valid())
			{
				std::stringstream ss; ss <<  " (Precond: " << m_spPrecond->name() << ")";
				convergence_check()->set_info(ss.str());
			}
			else
			{
				convergence_check()->set_info(" (No Preconditioner) ");
			}
		}

	///	writing debug output for a vector (if debug writer set)
		void write_debug_vector(const vector_type& vec, const char* filename)
		{
		//	if no debug writer set, we're done
			if(!m_pDebugWriter) return;

		//	write
			m_pDebugWriter->write_vector(vec, filename);
		}

	protected:
	// 	Operator that is inverted by this Inverse Operator
		ILinearOperator<vector_type,vector_type>* m_A;

	// 	Iterator used in the iterative scheme to compute the correction and update the defect
		SmartPtr<ILinearIterator<vector_type,vector_type> > m_spPrecond;

	//	flag if fresh defect should be computed when finish for debug purpose
		bool m_bRecomputeDefectWhenFinished;

	///	Debug Writer
		IDebugWriter<algebra_type>* m_pDebugWriter;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__ */

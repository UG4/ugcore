/*
 * linear_solver.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
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

template <typename TAlgebra>
class LinearSolver : public ILinearOperatorInverse<	typename TAlgebra::vector_type,
													typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	public:
		LinearSolver() :
			m_A(NULL), m_pPrecond(NULL), m_pConvCheck(NULL)
		{}

		LinearSolver( 	ILinearIterator<vector_type,vector_type>* Precond,
						IConvergenceCheck& ConvCheck) :
							m_A(NULL), m_pPrecond(Precond)
			{
				set_convergence_check(ConvCheck);
			};

		virtual const char* name() const {return "Iterative Linear Solver";}

		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
		}
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}
		void set_preconditioner(ILinearIterator<vector_type, vector_type>& precond)
		{
			m_pPrecond = &precond;
		}

		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
			m_A = &J;

			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(J, u))
				{
					UG_LOG("ERROR in 'LinearSolver::prepare': Cannot init "
							"Iterator Operator for Operator J.\n");return false;
				}

			return true;
		}

		virtual bool init(ILinearOperator<vector_type, vector_type>& L)
		{
			m_A = &L;

			// init Preconditioner for operator A
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(L))
				{
					UG_LOG("ERROR in 'LinearSolver::prepare': "
							"Cannot init Iterator Operator for Operator L.\n");
					return false;
				}

			return true;
		}

		virtual bool apply_return_defect(vector_type& cNLOut, vector_type& dNLIn)
		{
			if(m_A == NULL)
			{
				UG_LOG("ERROR: In 'LinearSolver::apply': Matrix A not set.\n");
				return false;
			}

			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'LinearSolver::apply': Convergence check not set.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!dNLIn.has_storage_type(PST_ADDITIVE) || !cNLOut.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("ERROR: In 'LinearSolver::apply':Inadequate storage format of Vectors.\n");
					return false;
				}
			#endif

			// rename d_nl as d (for convenience)
			vector_type& d = dNLIn;

			// build defect:  d := d_nl - J(u)*c_nl
			LS_PROFILE_BEGIN(LinSolve_BuildDefect);
			if(!m_A->apply_sub(d, cNLOut))
				{UG_LOG("ERROR in 'LinearSolver::apply': Unable to build defect. Aborting.\n"); return false;}
			LS_PROFILE_END(); //LinSolve_BuildDefect

			// create correction
			// todo: 	it would be sufficient to only copy the pattern (and parallel constructor)
			//			without initializing the values
			vector_type c; c.create(cNLOut.size()); c = cNLOut;

			prepare_conv_check();
			m_pConvCheck->start(d);

			// Iteration loop
			while(!m_pConvCheck->iteration_ended())
			{
				// Compute a correction c := B*d using one iterative step
				// Internally the defect is updated d := d - A*c = d - A*(x+c)
				if(m_pPrecond != NULL) {
					LS_PROFILE_BEGIN(LinSolve_PrecondApplyUpdateDefect);
					if(!m_pPrecond->apply_update_defect(c, d))
						{UG_LOG("ERROR in 'LinearSolver::apply': Iterator Operator "
									"applied incorrectly. Aborting.\n"); return false;}
					LS_PROFILE_END(); //LinSolve_PrecondApplyUpdateDefect
				}

				// add correction to solution
				LS_PROFILE_BEGIN(LinSolve_PrecondAddCorrection);
				cNLOut += c;
				LS_PROFILE_END(); //LinSolve_PrecondAddCorrection

				// compute new defect (in parallel)
				LS_PROFILE_BEGIN(LinSolve_ComputeNewDefect);
				m_pConvCheck->update(d);
				LS_PROFILE_END(); //LinSolve_ComputeNewDefect
			}

			LS_PROFILE_BEGIN(LinSolve_ConvCheckPost);
			if(!m_pConvCheck->post())
			{
				UG_LOG("ERROR in 'LinearSolver::apply': post-convergence-check signaled failure. Aborting.\n");
				return false;
			}
			LS_PROFILE_END(); //LinSolve_ConvCheckPost
			return true;
		}

		virtual bool apply(vector_type& cNLOut, const vector_type& dNLIn)
		{
		//	copy defect
			vector_type d; d.resize(dNLIn.size());
			d = dNLIn;

		//	solve on copy of defect
			return apply_return_defect(cNLOut, d);
		}

		// destructor
		virtual ~LinearSolver() {};

	protected:
		void prepare_conv_check()
		{
			m_pConvCheck->set_name(name());
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());
			if(m_pPrecond != NULL)
			{
				std::stringstream ss; ss <<  " (Precond: " << m_pPrecond->name() << ")";
				m_pConvCheck->set_info(ss.str());
			}
			else
			{
				m_pConvCheck->set_info(" (No Preconditioner) ");
			}
		}

	protected:
		// Operator that is inverted by this Inverse Operator
		ILinearOperator<vector_type,vector_type>* m_A;

		// Iterator used in the iterative scheme to compute the correction and update the defect
		ILinearIterator<vector_type,vector_type>* m_pPrecond;

		// Convergence Check
		IConvergenceCheck* m_pConvCheck;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__ */

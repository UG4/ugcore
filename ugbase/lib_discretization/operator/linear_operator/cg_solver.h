/*
 * cg_solver.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__

#include "lib_discretization/operator/operator.h"
#include "common/profiler/profiler.h"
#include "lib_discretization/io/vtkoutput.h"
namespace ug{

template <typename TFunction>
class CGSolver : public ILinearizedOperatorInverse<TFunction, TFunction>
{
	public:
		// domain space
		typedef TFunction domain_function_type;

		// range space
		typedef TFunction codomain_function_type;

	public:
			CGSolver( 	ILinearizedIteratorOperator<TFunction,TFunction>* Precond,
						ConvergenceCheck<TFunction>& ConvCheck) :
							m_pPrecond(Precond), m_ConvCheck(ConvCheck)
			{};

		virtual bool init(ILinearizedOperator<TFunction, TFunction>& A)
		{
			m_A = &A;

			// init Preconditioner for operator A
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(*m_A))
					{UG_LOG("ERROR in 'CGSolver::prepare': Cannot init Iterator Operator for Operator A.\n");return false;}

			return true;
		}

		// prepare Operator
		virtual bool prepare(TFunction& cNLOut, TFunction& uIn, TFunction& dNLIn)
		{
			m_pCurrentU = &uIn;

			return true;
		}

		// Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply(TFunction& xOut, TFunction& bIn)
		{
			#ifdef UG_PARALLEL
			if(!bIn.has_storage_type(PST_ADDITIVE) || !xOut.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("WARNING: In 'CGSolver::apply':Inadequate storage format of Vectors.\n");
					UG_LOG("                          use: b additive and x consistent to avoid internal type conversion.\n");
					if(!bIn.change_storage_type(PST_ADDITIVE)) return false;
					if(!xOut.change_storage_type(PST_CONSISTENT)) return false;
				}
			#endif

			// copy b as r
			domain_function_type& r = bIn;

			// build defect:  r := b - J(u)*x
			if(!m_A->apply_sub(r, xOut))
				{UG_LOG("ERROR in 'CGSolver::apply': Unable to build defect. Aborting.\n");return false;}

			// create help vector (h will be consistent r)
			domain_function_type t;
			codomain_function_type z, p;
			t.clone_pattern(r);
			z.clone_pattern(xOut);
			p.clone_pattern(xOut);

			// Preconditioning
			if(m_pPrecond != NULL)
			{
				if(!m_pPrecond->prepare(z, *m_pCurrentU, r))
					{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

				// apply z = M^-1 * s
				if(!m_pPrecond->apply(z, r, false))
					{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
			}
			else z = r;


			#ifdef UG_PARALLEL
			// make z consistent
			if(!z.change_storage_type(PST_CONSISTENT))
				{UG_LOG("Cannot convert z to consistent vector.\n"); return false;}
			#endif

			m_ConvCheck.set_offset(3);
			m_ConvCheck.set_symbol('%');
			m_ConvCheck.set_name("CG Solver");
			m_ConvCheck.start(r);

			number rho, rho_new, beta, alpha, lambda;
			rho = rho_new = beta = alpha = lambda = 0.0;

			// start search direction
			p = z;

			// start rho
			rho = VecProd(z, r);

			// Iteration loop
			while(!m_ConvCheck.iteration_ended())
			{
				// build t = A*p (t is additive afterwards)
				if(!m_A->apply(t, p))
					{UG_LOG("ERROR in 'CGSolver::apply': Unable "
								"to build t = A*p. Aborting.\n"); return false;}

				// compute alpha
				lambda = VecProd(t, p);
				alpha = rho/lambda;

				// update x := x + alpha*p
				//VecScaleAppend(xOut, p, alpha);
				x = x + alpha*p;


				// update r := r - alpha*t
				// VecScaleAppend(r, t, (-1)*alpha);
				r = r - alpha*t;

				// check convergence
				m_ConvCheck.update(r);

				// Preconditioning
				if(m_pPrecond != NULL)
				{
					if(!m_pPrecond->prepare(z, *m_pCurrentU, r))
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply z = M^-1 * r
					if(!m_pPrecond->apply(z, r, false))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
				}
				else z = r;


				#ifdef UG_PARALLEL
				// make z consistent
				if(!z.change_storage_type(PST_CONSISTENT))
					{UG_LOG("Cannot convert z to consistent vector.\n"); return false;}
				#endif

				// new rho
				rho_new = VecProd(z, r);

				// new beta
				beta = rho_new/rho;

				// new direction p:= beta*p + z
				p *= beta;
				p += z;

				// update rho
				rho = rho_new;
			}

			return m_ConvCheck.post();
		}

		// destructor
		virtual ~CGSolver() {};

	protected:
		number VecProd(domain_function_type& a, domain_function_type& b)
		{
			return a.dotprod(b);
		}

	protected:
		// Operator that is inverted by this Inverse Operator
		ILinearizedOperator<TFunction,TFunction>* m_A;

		// Iterator used in the iterative scheme to compute the correction and update the defect
		ILinearizedIteratorOperator<TFunction,TFunction>* m_pPrecond;

		// Convergence Check
		ConvergenceCheck<TFunction>& m_ConvCheck;

		// current solution
		TFunction* m_pCurrentU;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__ */

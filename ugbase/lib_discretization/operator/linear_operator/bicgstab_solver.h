/*
 * bicgstab_solver.h
 *
 *  Created on: 05.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__

#include "lib_discretization/operator/operator.h"
#include "common/profiler/profiler.h"

namespace ug{

template <typename TFunction>
class BiCGStabSolver : public ILinearizedOperatorInverse<TFunction, TFunction>
{
	public:
		// domain space
		typedef TFunction domain_function_type;

		// range space
		typedef TFunction codomain_function_type;

	public:
		BiCGStabSolver	( 	ILinearizedIteratorOperator<TFunction,TFunction>* Precond,
							ConvergenceCheck<TFunction>& ConvCheck) :
								m_pPrecond(Precond), m_ConvCheck(ConvCheck)
						{};

		virtual bool init(ILinearizedOperator<TFunction, TFunction>& A)
		{
			m_A = &A;

			// init Preconditioner for operator A
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(*m_A))
					{UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot init "
							"Iterator Operator for Operator A.\n"); return false;}

			return true;
		}

		// prepare Operator
		virtual bool prepare(TFunction& xOut, TFunction& uIn, TFunction& bIn)
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

			// build defect:  b := b - J(u)*x
			if(!m_A->apply_sub(bIn, xOut))
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
						"Unable to build defect. Aborting.\n");return false;}

			// create start r_0^* vector
			codomain_function_type r, p, v, q, t, s;
			r.clone_pattern(bIn);
			p.clone_pattern(bIn);
			v.clone_pattern(bIn);
			q.clone_pattern(xOut);
			t.clone_pattern(bIn);
			s.clone_pattern(bIn);

			m_ConvCheck.set_offset(3);
			m_ConvCheck.set_symbol('%');
			m_ConvCheck.set_name("BiCGStab Solver");
			m_ConvCheck.start(bIn);

			#ifdef UG_PARALLEL
			// convert b to unique (should already be unique due to norm calculation)
			if(!bIn.change_storage_type(PST_UNIQUE))
				{UG_LOG("Cannot convert b to unique vector.\n"); return false;}
			#endif

			number rho, rho_new, alpha, omega, beta, tt;
			tt = rho = alpha = omega = 0;

			// Iteration loop
			while(!m_ConvCheck.iteration_ended())
			{
				if(m_ConvCheck.step() == 0 /*or restart*/)
				{
					if(m_ConvCheck.step() != 0)
					{
						m_ConvCheck.update(bIn);
						if(m_ConvCheck.iteration_ended()) break;
					}

					// reset vectors
					r = bIn;

					// make r additive unique
				#ifdef UG_PARALLEL
					if(!r.change_storage_type(PST_UNIQUE))
						{UG_LOG("Cannot convert r to unique vector.\n"); return false;}
				#endif
					p.set(0.0);
					v.set(0.0);
					rho = alpha = omega = 1.0;
				}

				// compute rho new
				rho_new = VecProd(bIn, r);

				// compute new beta
				if(rho != 0.0 && omega != 0.0) beta = (rho_new/rho) * (alpha/omega);
				else {UG_LOG("rho= " << rho << " and omega= " << omega << " are invalid values. Aborting.\n"); return false;}

				// scale p := beta * p
				p *= beta;

				// add b to p (p:= p + b)
				p += bIn;

				// subtract: p := p - beta * omega * v
				VecScaleAppend(p, v, (-1)*beta*omega);


				// if preconditioner given
				if(m_pPrecond != NULL)
				{
					if(!m_pPrecond->prepare(q, *m_pCurrentU, p))
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply q = M^-1 * p
					if(!m_pPrecond->apply(q, p, false))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}

					// compute v := A*q
					if(!m_A->apply(v, q))
						{UG_LOG("ERROR: Unable to apply A. Aborting.\n");return false;}

					#ifdef UG_PARALLEL
					// make v unique
					if(!v.change_storage_type(PST_UNIQUE))
						{UG_LOG("Cannot convert v to unique vector.\n"); return false;}
					#endif

					alpha = VecProd(v, r);

					if(alpha != 0.0) alpha = rho_new/alpha;
					else {UG_LOG("alpha= " << alpha << " is an invalid value. Aborting.\n"); return false;}

					// add: x := x + alpha * q
					VecScaleAppend(xOut, q, alpha);
				}
				else
				{
					#ifdef UG_PARALLEL
					// make p consistent
					if(!p.change_storage_type(PST_CONSISTENT))
						{UG_LOG("Cannot convert p to consistent vector.\n"); return false;}
					#endif

					// compute v := A*p
					if(m_A->apply(v, p) != true)
						{UG_LOG("ERROR: Unable to apply A. Aborting.\n");return false;}

					#ifdef UG_PARALLEL
					// make v unique
					if(!v.change_storage_type(PST_UNIQUE))
						{UG_LOG("Cannot convert v to unique vector.\n"); return false;}
					#endif

					alpha = VecProd(v, r);

					if(alpha != 0.0) alpha = rho_new/alpha;
					else {UG_LOG("alpha= " << alpha << " is an invalid value. Aborting.\n"); return false;}

					// add: x := x + alpha * p
					VecScaleAppend(xOut, p, alpha);
				}


				// set s := b
				s = bIn;

				// update s := s - alpha*v
				VecScaleAppend(s, v, (-1)*alpha);

				// check convergence
				m_ConvCheck.update(s);
				if(m_ConvCheck.iteration_ended())
				{
					bIn = s;
					break;
				}

				// if preconditioner given
				if(m_pPrecond != NULL)
				{
					if(m_pPrecond->prepare(q, *m_pCurrentU, s) != true)
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply q = M^-1 * t
					if(!m_pPrecond->apply(q, s, false))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
				}
				else
				{
					// set q:=s
					q = s;

					#ifdef UG_PARALLEL
					// make q consistent
					if(!q.change_storage_type(PST_CONSISTENT))
						{UG_LOG("Cannot convert q to consistent vector.\n"); return false;}
					#endif
					}

				// compute t := A*q
				if(m_A->apply(t, q) != true)
					{UG_LOG("ERROR: Unable to apply A. Aborting.\n");return false;}

				#ifdef UG_PARALLEL
				// make t unique
				if(!t.change_storage_type(PST_UNIQUE))
					{UG_LOG("Cannot convert v to unique vector.\n"); return false;}
				#endif

				// tt = (t,t)
				tt = VecProd(t, t);

				// omega = (s,t)
				omega = VecProd(s, t);

				// omega = omega/tt
				if(tt != 0.0) omega = omega/tt;
				else {UG_LOG("tt= " << tt << " is an invalid value. Aborting.\n"); return false;}

				// add: x := x + omega * q
				VecScaleAppend(xOut, q, omega);

				// set b := s
				bIn = s;

				// 2. update of b:  b:= b - omega*t
				VecScaleAppend(bIn, t, (-1)*omega);

				// check convergence
				m_ConvCheck.update(bIn);

				// remember current rho
				rho = rho_new;
			}

			return m_ConvCheck.post();
		}

		// destructor
		virtual ~BiCGStabSolver() {};

	protected:
		bool VecScaleAppend(domain_function_type& a_func, domain_function_type& b_func, number s)
		{
			#ifdef UG_PARALLEL
			if(a_func.has_storage_type(PST_UNIQUE) && b_func.has_storage_type(PST_UNIQUE));
			else if(a_func.has_storage_type(PST_CONSISTENT) && b_func.has_storage_type(PST_CONSISTENT));
			else if (a_func.has_storage_type(PST_ADDITIVE) && b_func.has_storage_type(PST_ADDITIVE));
			else
			{
				a_func.change_storage_type(PST_ADDITIVE);
				b_func.change_storage_type(PST_ADDITIVE);
			}
			#endif
			typename domain_function_type::vector_type& a = a_func.get_vector();
			typename domain_function_type::vector_type& b = b_func.get_vector();

            for(size_t i = 0; i < a.size(); ++i)
				a[i] += s*b[i];
			return true;
		}

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

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__ */

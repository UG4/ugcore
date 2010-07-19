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
			return true;
		}

		// prepare Operator
		virtual bool prepare(codomain_function_type& u, domain_function_type& b, codomain_function_type& x)
		{
			// init iterator B for operator A
			if(m_pPrecond != NULL)
				if(m_pPrecond->init(*m_A) != true)
				{UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot init "
							"Iterator Operator for Operator A.\n"); return false;}

			m_pCurrentU = &u;

			return true;
		}

		// Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply(domain_function_type& b, codomain_function_type& x)
		{
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
			{UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
					"Wrong storage format of Vectors. Aborting.\n");return false;}

			// build defect:  b := b - J(u)*x
			if(m_A->apply_sub(x, b) != true)
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
						"Unable to build defect. Aborting.\n");return false;}

			// create start r_0^* vector
			codomain_function_type r, p, v, q, t, s;
			r.clone_pattern(b);
			p.clone_pattern(b);
			v.clone_pattern(b);
			q.clone_pattern(x);
			t.clone_pattern(b);
			s.clone_pattern(b);

			m_ConvCheck.set_offset(3);
			m_ConvCheck.set_symbol('%');
			m_ConvCheck.set_name("BiCGStab Solver");
			m_ConvCheck.start(b);

			// convert b to unique (should already be unique due to norm calculation)
			if(!b.change_storage_type(PST_UNIQUE))
				{UG_LOG("Cannot convert b to unique vector.\n"); return false;}

			number rho, rho_new, alpha, omega, beta, tt;
			tt = rho = alpha = omega = 0;

			// Iteration loop
			while(!m_ConvCheck.iteration_ended())
			{
				if(m_ConvCheck.step() == 0 /*or restart*/)
				{
					if(m_ConvCheck.step() != 0)
					{
						m_ConvCheck.update(b);
						if(m_ConvCheck.iteration_ended()) break;
					}

					// reset vectors
					r = b;

					// make r additive unique
					if(!r.change_storage_type(PST_UNIQUE))
						{UG_LOG("Cannot convert r to unique vector.\n"); return false;}
					p.set(0.0);
					v.set(0.0);
					rho = alpha = omega = 1.0;
				}

				// compute rho new
				rho_new = VecProd(b, r);

				// compute new beta
				if(rho != 0.0 && omega != 0.0) beta = (rho_new/rho) * (alpha/omega);
				else {UG_LOG("rho= " << rho << " and omega= " << omega << " are invalid values. Aborting.\n"); return false;}

				// scale p := beta * p
				p *= beta;

				// add b to p (p:= p + b)
				p += b;

				// subtract: p := p - beta * omega * v
				VecScaleAdd(p, v, (-1)*beta*omega);


				// if preconditioner given
				if(m_pPrecond != NULL)
				{
					// reset q
					q.set(0.0);

					// set s
					s = p;

					if(m_pPrecond->prepare(*m_pCurrentU, s, q) != true)
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply q = M^-1 * p
					if(!m_pPrecond->apply(s, q, true))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}

					// compute v := A*q
					if(!m_A->apply(q, v))
						{UG_LOG("ERROR: Unable to apply A. Aborting.\n");return false;}

					// make v unique
					if(!v.change_storage_type(PST_UNIQUE))
						{UG_LOG("Cannot convert v to unique vector.\n"); return false;}

					alpha = VecProd(v, r);

					if(alpha != 0.0) alpha = rho_new/alpha;
					else {UG_LOG("alpha= " << alpha << " is an invalid value. Aborting.\n"); return false;}

					// add: x := x + alpha * q
					VecScaleAdd(x, q, alpha);
				}
				else
				{
					// make p consistent
					if(!p.change_storage_type(PST_CONSISTENT))
						{UG_LOG("Cannot convert p to consistent vector.\n"); return false;}

					// compute v := A*p
					if(m_A->apply(p, v) != true)
						{UG_LOG("ERROR: Unable to apply A. Aborting.\n");return false;}

					// make v unique
					if(!v.change_storage_type(PST_UNIQUE))
						{UG_LOG("Cannot convert v to unique vector.\n"); return false;}

					alpha = VecProd(v, r);

					if(alpha != 0.0) alpha = rho_new/alpha;
					else {UG_LOG("alpha= " << alpha << " is an invalid value. Aborting.\n"); return false;}

					// add: x := x + alpha * p
					VecScaleAdd(x, p, alpha);
				}


				// set s := b
				s = b;

				// update s := s - alpha*v
				VecScaleAdd(s, v, (-1)*alpha);

				// check convergence
				m_ConvCheck.update(s);
				if(m_ConvCheck.iteration_ended())
				{
					b = s;
					break;
				}

				// if preconditioner given
				if(m_pPrecond != NULL)
				{
					// reset q
					q.set(0.0);

					// copy s
					t = s;

					if(m_pPrecond->prepare(*m_pCurrentU, t, q) != true)
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply q = M^-1 * t
					if(!m_pPrecond->apply(t, q, true))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
				}
				else
				{
					// set q:=s
					q = s;

					// make q consistent
					if(!q.change_storage_type(PST_CONSISTENT))
						{UG_LOG("Cannot convert q to consistent vector.\n"); return false;}
				}

				// compute t := A*q
				if(m_A->apply(q, t) != true)
					{UG_LOG("ERROR: Unable to apply A. Aborting.\n");return false;}

				// make t unique
				if(!t.change_storage_type(PST_UNIQUE))
					{UG_LOG("Cannot convert v to unique vector.\n"); return false;}

				// tt = (t,t)
				tt = VecProd(t, t);

				// omega = (s,t)
				omega = VecProd(s, t);

				// omega = omega/tt
				if(tt != 0.0) omega = omega/tt;
				else {UG_LOG("tt= " << tt << " is an invalid value. Aborting.\n"); return false;}

				// add: x := x + omega * q
				VecScaleAdd(x, q, omega);

				// set b := s
				b = s;

				// 2. update of b:  b:= b - omega*t
				VecScaleAdd(b, t, (-1)*omega);

				// check convergence
				m_ConvCheck.update(b);

				// remember current rho
				rho = rho_new;
			}

			return m_ConvCheck.post();
		}

		// destructor
		virtual ~BiCGStabSolver() {};

	protected:
		bool VecScaleAdd(domain_function_type& a_func, domain_function_type& b_func, number s)
		{
			typename domain_function_type::vector_type& a = a_func.get_vector();
			typename domain_function_type::vector_type& b = b_func.get_vector();
			typename domain_function_type::algebra_type::matrix_type::local_matrix_type locMat(1, 1);
			typename domain_function_type::algebra_type::matrix_type::local_index_type locInd(1);
			typename domain_function_type::vector_type::local_vector_type locVec(1);

			for(size_t i = 0; i < a.size(); ++i){
				locInd[0][0] = i;
				b.get(locVec, locInd);
				locVec[0] *= s;

				a.add(locVec, locInd);
			}
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

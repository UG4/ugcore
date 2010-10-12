/*
 * bicgstab_solver.h
 *
 *  Created on: 05.07.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__

#include "lib_algebra/operator/operator_interface.h"
#include "common/profiler/profiler.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TAlgebra>
class BiCGStabSolver : public ILinearOperatorInverse< 	typename TAlgebra::vector_type,
														typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	public:
		BiCGStabSolver() :
			m_pPrecond(NULL), m_pConvCheck(NULL)
		{};

		BiCGStabSolver	( 	ILinearIterator<vector_type,vector_type>* Precond,
							IConvergenceCheck& ConvCheck) :
							m_pPrecond(Precond), m_pConvCheck(&ConvCheck)
						{};

		void set_convergence_check(IConvergenceCheck& convCheck) {m_pConvCheck = &convCheck;}
		void set_preconditioner(ILinearIterator<vector_type, vector_type>& precond) {m_pPrecond = &precond;}

		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
			m_A = &J;

			// init Preconditioner for operator J
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(J, u))
				{
					UG_LOG("ERROR in 'BiCGStabSolver::prepare': Cannot init "
							"Iterator Operator for Operator J.\n"); return false;
				}

			return true;
		}

		virtual bool init(ILinearOperator<vector_type, vector_type>& L)
		{
			m_A = &L;

			// init Preconditioner for operator L
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(L))
				{
					UG_LOG("ERROR in 'CGSolver::prepare': "
							"Cannot init Iterator Operator for Operator L.\n");
					return false;
				}

			return true;
		}

		// Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& xOut, vector_type& bIn)
		{
			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'BiCGStabSolver::apply': Convergence check not set.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!bIn.has_storage_type(PST_ADDITIVE) || !xOut.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("ERROR: In 'BiCGStabSolver::apply':Inadequate storage format of Vectors.\n");
					return false;
				}
			#endif

			// build defect:  b := b - J(u)*x
			if(!m_A->apply_sub(bIn, xOut))
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': "
						"Unable to build defect. Aborting.\n");return false;}

			// create start r_0^* vector
		//	todo: 	It would be sufficient to copy only the pattern and
		//			without initializing
			vector_type r; r.create(bIn.size()); r = bIn;
			vector_type p; p.create(bIn.size()); p = bIn;
			vector_type v; v.create(bIn.size()); v = bIn;
			vector_type t; t.create(bIn.size()); t = bIn;
			vector_type s; s.create(bIn.size()); s = bIn;
			vector_type q; q.create(xOut.size()); q = xOut;

			m_pConvCheck->set_offset(3);
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name("BiCGStab Solver");
			m_pConvCheck->start(bIn);

			#ifdef UG_PARALLEL
			// convert b to unique (should already be unique due to norm calculation)
			if(!bIn.change_storage_type(PST_UNIQUE))
				{UG_LOG("Cannot convert b to unique vector.\n"); return false;}
			#endif

			number rho, rho_new, alpha, omega, beta, tt;
			tt = rho = alpha = omega = 0;

		// 	Iteration loop
			while(!m_pConvCheck->iteration_ended())
			{
				if(m_pConvCheck->step() == 0 /*or restart*/)
				{
					if(m_pConvCheck->step() != 0)
					{
						m_pConvCheck->update(bIn);
						if(m_pConvCheck->iteration_ended()) break;
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

			// 	Compute rho new
				rho_new = VecProd(bIn, r);

			// 	Compute new beta
				if(rho != 0.0 && omega != 0.0) beta = (rho_new/rho) * (alpha/omega);
				else {UG_LOG("rho= " << rho << " and omega= " << omega << " are invalid values. Aborting.\n"); return false;}

			// 	scale p := beta * p
				p *= beta;

			// 	add b to p (p:= p + b)
				p += bIn;

			// 	subtract: p := p - beta * omega * v
				VecScaleAppend(p, v, (-1)*beta*omega);


			// 	if preconditioner given
				if(m_pPrecond != NULL)
				{
				// 	apply q = M^-1 * p
					if(!m_pPrecond->apply(q, p))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}

				// 	compute v := A*q
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

				// 	add: x := x + alpha * q
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
				m_pConvCheck->update(s);
				if(m_pConvCheck->iteration_ended())
				{
					bIn = s;
					break;
				}

				// if preconditioner given
				if(m_pPrecond != NULL)
				{
					// apply q = M^-1 * t
					if(!m_pPrecond->apply(q, s))
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
				m_pConvCheck->update(bIn);

				// remember current rho
				rho = rho_new;
			}

			return m_pConvCheck->post();
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
		virtual ~BiCGStabSolver() {};

	protected:
		bool VecScaleAppend(vector_type& a, vector_type& b, number s)
		{
			#ifdef UG_PARALLEL
			if(a.has_storage_type(PST_UNIQUE) && b.has_storage_type(PST_UNIQUE));
			else if(a.has_storage_type(PST_CONSISTENT) && b.has_storage_type(PST_CONSISTENT));
			else if (a.has_storage_type(PST_ADDITIVE) && b.has_storage_type(PST_ADDITIVE));
			else
			{
				a.change_storage_type(PST_ADDITIVE);
				b.change_storage_type(PST_ADDITIVE);
			}
			#endif

            for(size_t i = 0; i < a.size(); ++i)
				a[i] += s*b[i];
			return true;
		}

		number VecProd(vector_type& a, vector_type& b)
		{
			return a.dotprod(b);
		}

	protected:
		// Operator that is inverted by this Inverse Operator
		ILinearOperator<vector_type,vector_type>* m_A;

		// Iterator used in the iterative scheme to compute the correction and update the defect
		ILinearIterator<vector_type,vector_type>* m_pPrecond;

		// Convergence Check
		IConvergenceCheck* m_pConvCheck;

		// current solution
		vector_type* m_pCurrentU;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__ */

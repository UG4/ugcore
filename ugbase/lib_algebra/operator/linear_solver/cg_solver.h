/*
 * cg_solver.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__
#include <iostream>
#include <sstream>

#include "lib_algebra/operator/operator_interface.h"
#include "common/profiler/profiler.h"
#ifdef UG_PARALLEL
	#include "lib_algebra/parallelization/parallelization.h"
#endif

namespace ug{

template <typename TAlgebra>
class CGSolver : public ILinearOperatorInverse<	typename TAlgebra::vector_type,
													typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	public:
		CGSolver() :
			m_pPrecond(NULL), m_pConvCheck(NULL)
		{}

		CGSolver( 	ILinearIterator<vector_type,vector_type>* Precond,
						IConvergenceCheck& ConvCheck) :
							m_pPrecond(Precond), m_pConvCheck(&ConvCheck)
			{};

		virtual const char* name() const {return "BiCGStabSolver";}

		void set_convergence_check(IConvergenceCheck& convCheck)
		{
			m_pConvCheck = &convCheck;
			m_pConvCheck->set_offset(3);
			m_pConvCheck->set_symbol('%');
			m_pConvCheck->set_name(name());

			if(m_pPrecond != NULL)
			{
				stringstream ss; ss <<  " (Precond: " << m_pPrecond->name() << ")";
				m_pConvCheck->set_info(ss.str());
			}
		}
		IConvergenceCheck* get_convergence_check() {return m_pConvCheck;}
		void set_preconditioner(ILinearIterator<vector_type, vector_type>& precond)
		{
			m_pPrecond = &precond;
			if(m_pConvCheck != NULL)
			{
				stringstream ss; ss <<  " (Precond: " << m_pPrecond->name() << ")";
				m_pConvCheck->set_info(ss.str());
			}
		}

		virtual bool init(ILinearOperator<vector_type, vector_type>& J, const vector_type& u)
		{
			m_A = &J;

			// init Preconditioner for operator A
			if(m_pPrecond != NULL)
				if(!m_pPrecond->init(J, u))
				{
					UG_LOG("ERROR in 'CGSolver::prepare': "
							"Cannot init Iterator Operator for Operator A.\n");
					return false;
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
					UG_LOG("ERROR in 'CGSolver::prepare': "
							"Cannot init Iterator Operator for Operator A.\n");
					return false;
				}

			return true;
		}

		// Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply_return_defect(vector_type& xOut, vector_type& bIn)
		{
			if(m_pConvCheck == NULL)
			{
				UG_LOG("ERROR: In 'CGSolver::apply': Convergence check not set.\n");
				return false;
			}

			#ifdef UG_PARALLEL
			if(!bIn.has_storage_type(PST_ADDITIVE) || !xOut.has_storage_type(PST_CONSISTENT))
				{
					UG_LOG("ERROR: In 'CGSolver::apply':Inadequate storage format of Vectors.\n");
					return false;
				}
			#endif

		// 	rename r as b (for convenience)
			vector_type& r = bIn;

		// 	Build defect:  r := b - J(u)*x
			if(!m_A->apply_sub(r, xOut))
				{UG_LOG("ERROR in 'CGSolver::apply': "
						"Unable to build defect. Aborting.\n");return false;}

		// 	create help vector (h will be consistent r)
		//	todo: 	It would be sufficient to copy only the pattern and
		//			without initializing, but in parallel we have to copy communicators
			vector_type t; t.create(r.size()); t = r;
			vector_type z; z.create(xOut.size()); z = xOut;
			vector_type p; p.create(xOut.size()); p = xOut;

		// 	Preconditioning
			if(m_pPrecond != NULL)
			{
				// apply z = M^-1 * s
				if(!m_pPrecond->apply(z, r))
					{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
			}
			else z = r;


			#ifdef UG_PARALLEL
			// make z consistent
			if(!z.change_storage_type(PST_CONSISTENT))
				{UG_LOG("Cannot convert z to consistent vector.\n"); return false;}
			#endif

			m_pConvCheck->start(r);

			number rho, rho_new, beta, alpha, lambda;
			rho = rho_new = beta = alpha = lambda = 0.0;

			// start search direction
			p = z;

			// start rho
			rho = VecProd(z, r);

		// 	Iteration loop
			while(!m_pConvCheck->iteration_ended())
			{
			// 	Build t = A*p (t is additive afterwards)
				if(!m_A->apply(t, p))
					{UG_LOG("ERROR in 'CGSolver::apply': Unable "
								"to build t = A*p. Aborting.\n"); return false;}

			// 	Compute alpha
				lambda = VecProd(t, p);
				alpha = rho/lambda;

			// 	Update xOut := xOut + alpha*p
				VecScaleAdd(xOut, 1.0, xOut, alpha, p);

			// 	Update r := r - alpha*t
				VecScaleAdd(r, 1.0, r, -alpha, t);

			// 	Check convergence
				m_pConvCheck->update(r);

			// 	Preconditioning
				if(m_pPrecond != NULL)
				{
				// 	apply z = M^-1 * r
					if(!m_pPrecond->apply(z, r))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
				}
				else z = r;


				#ifdef UG_PARALLEL
			// 	make z consistent
				if(!z.change_storage_type(PST_CONSISTENT))
					{UG_LOG("Cannot convert z to consistent vector.\n"); return false;}
				#endif

			// 	new rho
				rho_new = VecProd(z, r);

			// 	new beta
				beta = rho_new/rho;

			// 	new direction p:= beta*p + z
				VecScaleAdd(p, beta, p, 1.0, z);

			// 	update rho
				rho = rho_new;
			}

			return m_pConvCheck->post();
		}

		virtual bool apply(vector_type& cNLOut, const vector_type& dNLIn)
		{
		//	copy defect
			vector_type d;
			d = dNLIn;

		//	solve on copy of defect
			return apply_return_defect(cNLOut, d);
		}

		// destructor
		virtual ~CGSolver() {};

	protected:
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

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__ */

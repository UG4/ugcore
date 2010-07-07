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

namespace ug{

template <typename TDiscreteFunction>
class CGSolver : public ILinearizedOperatorInverse<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// domain space
		typedef TDiscreteFunction domain_function_type;

		// range space
		typedef TDiscreteFunction codomain_function_type;

	public:
		//TODO: Implement usage of preconditioner
		CGSolver( 	ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* Precond,
					int maxIter, number absTol, number relTol,
					int verboseLevel = 0) :
			m_verboseLevel(verboseLevel), m_pIter(Precond),
			m_maxIter(maxIter), m_absTol(absTol), m_relTol(relTol)
			{};

		virtual bool init(ILinearizedOperator<TDiscreteFunction, TDiscreteFunction>& A)
		{
			m_A = &A;
			return true;
		}

		// prepare Operator
		virtual bool prepare(codomain_function_type& u, domain_function_type& d_nl, codomain_function_type& c_nl)
		{
			// init iterator B for operator A
			if(m_pIter != NULL)
				if(m_pIter->init(*m_A) != true)
				{UG_LOG("ERROR in 'CGSolver::prepare': Cannot init Iterator Operator for Operator A.\n");return false;}

			// prepare iterator B for d_nl and c_nl
			if(m_pIter != NULL)
				if(m_pIter->prepare(u, d_nl, c_nl) != true)
				{UG_LOG("ERROR in 'CGSolver::prepare': Cannot prepare Iterator Operator.\n"); return false;}

			m_pCurrentU = &u;

			return true;
		}

		// Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply(domain_function_type& b, codomain_function_type& x)
		{
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
				{UG_LOG("ERROR in 'CGSolver::apply': Wrong storage format of Vectors. Aborting.\n"); return false;}

			// copy b as r
			domain_function_type& r = b;

			// build defect:  r := b - J(u)*x
			if(m_A->apply_sub(x, r) != true)
				{UG_LOG("ERROR in 'CGSolver::apply': Unable to build defect. Aborting.\n");return false;}

			// create help vector (h will be consistent r)
			domain_function_type t;
			codomain_function_type z, p;
			t.clone_pattern(r);
			z.clone_pattern(x);
			p.clone_pattern(x);

			// Preconditioning
			if(m_pIter != NULL)
			{
				// copy r
				t = r;

				if(m_pIter->prepare(*m_pCurrentU, t, z) != true)
					{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

				// apply z = M^-1 * s
				if(!m_pIter->apply(t, z))
					{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
			}
			else
			{
				z = r;

				// make z consistent
				if(!z.change_storage_type(PST_CONSISTENT))
					{UG_LOG("Cannot convert z to consistent vector.\n"); return false;}
			}

			// compute start norm ||r||_2
			number norm, norm_old, norm_start;

			// compute first norm
			PROFILE_BEGIN(GlobalVecProd_FirstNorm);
				norm = norm_old = norm_start = r.two_norm();
			PROFILE_END();

			// Print Start information
			if(m_verboseLevel >= 1) UG_LOG("\n    %%%%%%%%%% CG Solver %%%%%%%%%%\n");
			if(m_verboseLevel >= 2) UG_LOG("    %   Iter     Defect         Rate \n");
			if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << 0 << ":  " << std::scientific << norm_old <<  "      -------" << std::endl);

			// exit flag
			bool bExit;

			// check convergence
			if(!check_convergence(norm, norm_start, 0, bExit)) return bExit;

			number rho, rho_new, beta, alpha, lambda;
			rho = rho_new = beta = alpha = lambda = 0.0;

			// start rho
			rho = VecProd(z, r);

			// start search direction
			p = z;

			// Iteration loop
			for(int i = 1; ; ++i)
			{
				// build t = A*p (t is additive afterwards)
				if(m_A->apply(p, t) != true)
					{UG_LOG("ERROR in 'CGSolver::apply': Unable to build t = A*p. Aborting.\n");return false;}

				// compute alpha
				lambda = VecProd(t, p);
				alpha = rho/lambda;

				// update x := x + alpha*p
				VecScaleAdd(x, p, alpha);

				// update r := r - alpha*t
				VecScaleAdd(r, t, (-1)*alpha);

				// compute new norm
				norm = r.two_norm();

				// print convergence rate
				if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// check convergence
				if(!check_convergence(norm, norm_start, i, bExit)) return bExit;

				// remember current norm
				norm_old = norm;

				// Preconditioning
				if(m_pIter != NULL)
				{
					// copy r
					t = r;

					if(m_pIter->prepare(*m_pCurrentU, t, z) != true)
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply z = M^-1 * s
					if(!m_pIter->apply(t, z))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}
				}
				else
				{
					z = r;

					// make z consistent
					if(!z.change_storage_type(PST_CONSISTENT))
						{UG_LOG("Cannot convert z to consistent vector.\n"); return false;}
				}

				// new rho
				rho_new = VecProd(z, r);

				// new beta
				beta = rho_new/rho;

				// new direction p:= beta*p + z
				p *= beta;
				p+= z;

				// update rho
				rho = rho_new;
			}
			UG_ASSERT(0, "This line should never be reached.");
			return false;
		}

		// destructor
		virtual ~CGSolver() {};

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
		bool check_convergence(number norm, number norm_start, int num_iter, bool& bExit)
		{
			// check that defect is a still a valid number
			if(!is_valid_number(norm))
			{
				if(m_verboseLevel >= 1) UG_LOG("    %%%%% Defect " << norm << " is not a valid number. Linear Solver did NOT CONVERGE. %%%%%\n\n");
				bExit = false;
				return false;
			}

			// check if defect is small enough (absolute)
			if(norm < m_absTol)
			{
				if(m_verboseLevel >= 1) UG_LOG("    %%%%% Absolute defect " << m_absTol << " reached. Linear Solver converged. %%%%%\n\n");
				bExit = true;
				return false;
			}

			// check if defect is small enough (relative)
			if(norm/norm_start < m_relTol)
			{
				if(m_verboseLevel >= 1) UG_LOG("    %%%%% Relative defect " << m_relTol << " reached. Linear Solver converged. %%%%%\n\n");
				bExit = true;
				return false;
			}

			// check that maximum number of iterations is not reached
			if(num_iter > m_maxIter)
			{
				if(m_verboseLevel >= 1) UG_LOG("    %%%%% Absolute defect " << m_absTol << " and relative defect " << m_relTol << " NOT reached after " << m_maxIter << " Iterations. %%%%%\n");
				if(m_verboseLevel >= 1) UG_LOG("    %%%%% Iterative Linear Solver did NOT CONVERGE. %%%%%\n\n");
				bExit = false;
				return false;
			}

			// return false, if still continuing iteration
			return true;
		}

		void print(int verboseLevel, std::ostream outStream)
		{
			if(verboseLevel >= m_verboseLevel) UG_LOG(outStream);
		}

		bool is_valid_number(number value)
		{
			// (value >= std::numeric_limits<number>::min() ) == true if value > -infty
			// (value <= std::numeric_limits<number>::max() ) == true if value < infty
			// (value == value                         ) == true if value != NaN

			if (value == 0.0) return true;
			else return value >= std::numeric_limits<number>::min() && value <= std::numeric_limits<number>::max() && value == value && value >= 0.0;
		}

		// Discribes, how many output is printed. (0 = nothing, 1 = major informations, 2 = all)
		int m_verboseLevel;

	protected:
		// Operator that is inverted by this Inverse Operator
		ILinearizedOperator<TDiscreteFunction,TDiscreteFunction>* m_A;

		// Iterator used in the iterative scheme to compute the correction and update the defect
		ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* m_pIter;

		TDiscreteFunction* m_pCurrentU;

		// maximal number of iterations
		int m_maxIter;

		// absolute defect to be reached
		number m_absTol;

		// relative defect to be reached
		number m_relTol;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__ */

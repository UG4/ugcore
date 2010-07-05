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

template <typename TDiscreteFunction>
class BiCGStabSolver : public ILinearizedOperatorInverse<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// domain space
		typedef TDiscreteFunction domain_function_type;

		// range space
		typedef TDiscreteFunction codomain_function_type;

	public:
		BiCGStabSolver( 	ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* Precond,
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
		virtual bool prepare(codomain_function_type& u, domain_function_type& b, codomain_function_type& x)
		{
			// TODO: Do we assume, that operator has been prepared? Do we have to prepare it here?
			// m_A->prepare(u, d_nl, c_nl);

			// init iterator B for operator A
			if(m_pIter != NULL)
				if(m_pIter->init(*m_A) != true)
				{UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot init Iterator Operator for Operator A.\n"); return false;}

			m_pCurrentU = &u;

			return true;
		}

		// Solve J(u)*x = b, such that x = J(u)^{-1} b
		virtual bool apply(domain_function_type& b, codomain_function_type& x)
		{
			if(!b.has_storage_type(PST_ADDITIVE) || !x.has_storage_type(PST_CONSISTENT))
			{
				UG_LOG("ERROR in 'LinearOperatorInverse::apply': Wrong storage format of Vectors. Aborting.\n");
				return false;
			}

			// build defect:  b := b - J(u)*x
			if(m_A->apply_sub(x, b) != true)
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to build defect. Aborting.\n");return false;}

			// create start r_0^* vector
			codomain_function_type r, p, v, q, t, s;
			r.clone_pattern(b);
			p.clone_pattern(b);
			v.clone_pattern(b);
			q.clone_pattern(x);
			t.clone_pattern(b);
			s.clone_pattern(b);

			// compute start norm ||d||_2
			number norm, norm_old, norm_start;

			// compute first norm
			PROFILE_BEGIN(GlobalVecProd_FirstNorm);
				norm = norm_old = norm_start = b.two_norm();
			PROFILE_END();

			// convert b to unique (should already be unique due to norm calculation)
			if(!b.change_storage_type(PST_UNIQUE))
				{UG_LOG("Cannot convert b to unique vector.\n"); return false;}

			// Print Start information
			if(m_verboseLevel >= 1) UG_LOG("\n    %%%%%%%%%% BiCGStab Solver %%%%%%%%%%\n");
			if(m_verboseLevel >= 2) UG_LOG("    %   Iter     Defect         Rate \n");
			if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << 0 << ":  " << std::scientific << norm_old <<  "      -------" << std::endl);

			number rho, rho_new, alpha, omega, beta, tt;
			tt = rho = alpha = omega = 0;
			bool bExit;

			// Iteration loop
			for(int i = 1; ; ++i)
			{
				if(i==1 /*or restart*/)
				{
					// check if already ready before iterating
					if(!check_convergence(norm, norm_start, i, bExit)) return bExit;

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
				if(m_pIter != NULL)
				{
					// reset q
					q.set(0.0);

					// set s
					s = p;

					if(m_pIter->prepare(*m_pCurrentU, s, q) != true)
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply q = M^-1 * p
					if(!m_pIter->apply(s, q))
						{UG_LOG("ERROR: Cannot apply preconditioner. Aborting.\n"); return false;}

					// compute v := A*q
					if(m_A->apply(q, v) != true)
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

					// add: x := x + alpha * q
					VecScaleAdd(x, p, alpha);
				}


				// set s := b
				s = b;

				// update s := s - alpha*v
				VecScaleAdd(s, v, (-1)*alpha);

				// compute norm
				norm = s.two_norm();

				// print convergence rate
				if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// remember current norm
				norm_old = norm; i++;

				// check if already ready before iterating
				if(!check_convergence(norm, norm_start, i, bExit))
				{
					b = s;
					return bExit;
				}


				// if preconditioner given
				if(m_pIter != NULL)
				{
					// reset q
					q.set(0.0);

					// copy s
					t = s;

					if(m_pIter->prepare(*m_pCurrentU, t, q) != true)
						{UG_LOG("ERROR: Cannot prepare preconditioner. Aborting.\n"); return false;}

					// apply q = M^-1 * p
					if(!m_pIter->apply(t, q))
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
				VecScaleAdd(x, q, alpha);

				// set b := s
				b = s;

				// 2. update of b:  b:= b - omega*t
				VecScaleAdd(b, t, (-1)*omega);

				// compute norm
				norm = b.two_norm();

				// print convergence rate
				if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// remember current norm
				norm_old = norm;

				// check if already ready before iterating
				if(!check_convergence(norm, norm_start, i, bExit)) return bExit;

				// remember current rho
				rho = rho_new;
			}
			UG_ASSERT(0, "This line should never be reached.");
			return false;
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
			bool check = false;
			if(a.has_storage_type(PST_ADDITIVE) && b.has_storage_type(PST_CONSISTENT)) check = true;
			if(b.has_storage_type(PST_ADDITIVE) && a.has_storage_type(PST_CONSISTENT)) check = true;
			if(b.has_storage_type(PST_UNIQUE) && a.has_storage_type(PST_UNIQUE)) check = true;

			if(!check)
			{
				// fall back, should be improved. For this function, we always expect to match the upper 3 cases
				a.change_storage_type(PST_ADDITIVE);
				b.change_storage_type(PST_CONSISTENT);
			}

			typename domain_function_type::vector_type& a_vec = a.get_vector();
			typename domain_function_type::vector_type& b_vec = b.get_vector();

			// step 1: compute local dot product
			double tSumLocal = (double)a_vec.dotprod(b_vec);
			double tSumGlobal;

			// step 2: compute new defect norm
#ifdef UG_PARALLEL
			pcl::AllReduce(&tSumLocal, &tSumGlobal, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
#else
			tSumGlobal = tSumLocal;
#endif

			return tSumGlobal;
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
				if(m_verboseLevel >= 1) UG_LOG("    %%%%% Absolute defect " << m_absTol << " and relative defect " << m_relTol << " NOT reached after " << m_maxIter << " Iterations. %%%%%");
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

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__BICGSTAB_SOLVER__ */

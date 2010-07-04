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
		CGSolver( 	ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>& B,
					int maxIter, number absTol, number relTol,
					int verboseLevel = 0) :
			m_verboseLevel(verboseLevel), m_iter(&B),
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
			// TODO: Do we assume, that operator has been prepared? Do we have to prepare it here?
			// m_A->prepare(u, d_nl, c_nl);

			// init iterator B for operator A
			if(m_iter->init(*m_A) != true)
			{
				UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot init Iterator Operator for Operator A.\n");
				return false;
			}

			// prepare iterator B for d_nl and c_nl
			if(m_iter->prepare(u, d_nl, c_nl) != true)
			{
				UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot prepare Iterator Operator.\n");
				return false;
			}

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

			// copy b as r
			domain_function_type& r = b;

			// build defect:  r := b - J(u)*x
			if(m_A->apply_sub(x, r) != true)
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to build defect. Aborting.\n");return false;}

			// create help vector (h will be consistent r)
			codomain_function_type h(r);

			// make h consistent
			if(!h.change_storage_type(PST_CONSISTENT))
				{UG_LOG("Cannot convert h to consistent vector.\n"); return false;}

			// create help vector (search direction)
			codomain_function_type d(h);

			// compute start norm ||d||_2
			number norm, norm_old, norm_start;

			// compute first norm
			PROFILE_BEGIN(GlobalVecProd_FirstNorm);
				norm = norm_old = norm_start = VecProd(h, r);
			PROFILE_END();

			// Print Start information
			if(m_verboseLevel >= 1) UG_LOG("\n    %%%%%%%%%% CG Solver %%%%%%%%%%\n");
			if(m_verboseLevel >= 2) UG_LOG("    %   Iter     Defect         Rate \n");
			if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << 0 << ":  " << std::scientific << norm_old <<  "      -------" << std::endl);

			// Iteration loop
			for(int i = 1; ; ++i)
			{
				// check that defect is a still a valid number
				if(!is_valid_number(norm))
				{
					if(m_verboseLevel >= 1) UG_LOG("    %%%%% Defect " << norm << " is not a valid number. Linear Solver did NOT CONVERGE. %%%%%\n\n");
					return false;
				}

				// check if defect is small enough (absolute)
				if(norm < m_absTol)
				{
					if(m_verboseLevel >= 1) UG_LOG("    %%%%% Absolute defect " << m_absTol << " reached. Linear Solver converged. %%%%%\n\n");
					return true;
				}

				// check if defect is small enough (relative)
				if(norm/norm_start < m_relTol)
				{
					if(m_verboseLevel >= 1) UG_LOG("    %%%%% Relative defect " << m_relTol << " reached. Linear Solver converged. %%%%%\n\n");
					return true;
				}

				// check that maximum number of iterations is not reached
				if(i > m_maxIter)
				{
					if(m_verboseLevel >= 1) UG_LOG("    %%%%% Absolute defect " << m_absTol << " and relative defect " << m_relTol << " NOT reached after " << m_maxIter << " Iterations. %%%%%");
					if(m_verboseLevel >= 1) UG_LOG("    %%%%% Iterative Linear Solver did NOT CONVERGE. %%%%%\n\n");
					return false;
				}

				// build h = A*d (h is additive afterwards)
				if(m_A->apply(d, h) != true)
					{UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to build h = A*d. Aborting.\n");return false;}

				// compute alpha
				number alpha = VecProd(d, h);
				alpha = norm/alpha;

				// update x := x + alpha*d
				VecScaleAdd(x, d, alpha);

				// update r := r - alpha*h
				VecScaleAdd(r, h, (-1)*alpha);

				// copy new r into h
				h = r;

				// make h consistent
				if(!h.change_storage_type(PST_CONSISTENT))
					{UG_LOG("Cannot convert h to consistent vector.\n"); return false;}

				// compute new norm
				PROFILE_BEGIN(GlobalVecProd_NewNorm);
				norm = VecProd(h, r);
				PROFILE_END();

				// new beta
				number beta = norm/norm_old;

				// new direction d:= beta*d + h
				d *= beta;
				d += h;

				// print convergence rate
				if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// remember current norm
				norm_old = norm;
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
			bool check = false;
			if(a.has_storage_type(PST_ADDITIVE) && b.has_storage_type(PST_CONSISTENT)) check = true;
			if(b.has_storage_type(PST_ADDITIVE) && a.has_storage_type(PST_CONSISTENT)) check = true;

			if(!check) return -1;

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
		ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>* m_iter;

		// maximal number of iterations
		int m_maxIter;

		// absolute defect to be reached
		number m_absTol;

		// relative defect to be reached
		number m_relTol;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__CG_SOLVER__ */

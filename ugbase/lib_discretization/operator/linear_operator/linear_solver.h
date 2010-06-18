/*
 * operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__

#include "lib_discretization/operator/operator.h"

namespace ug{

template <typename TDiscreteFunction>
class LinearSolver : public ILinearizedOperatorInverse<TDiscreteFunction, TDiscreteFunction>
{
	public:
		// domain space
		typedef TDiscreteFunction domain_function_type;

		// range space
		typedef TDiscreteFunction codomain_function_type;

	public:
		LinearSolver( 	ILinearizedIteratorOperator<TDiscreteFunction,TDiscreteFunction>& B,
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

		// Solve J(u)*c_nl = d_nl, such that c_nl = J(u)^{-1} d_nl
		// This is done by iterating: c_nl := c_nl + B(u)(d_nl - J(u)*c_nl)
		// In d_nl the last defect d := d_nl - J(u)*c_nl is returned
		// In the following:
		// c_nl, d_nl refer to the non-linear defect and correction as e.g. in J(u) * c_nl = d_nl as it appears in Newton scheme
		// c, d are the correction and defect for solving that linear equation iteratively.
		virtual bool apply(domain_function_type& d_nl, codomain_function_type& c_nl)
		{
			if(!d_nl.has_storage_type(GFST_ADDITIVE) || !c_nl.has_storage_type(GFST_CONSISTENT))
			{
				UG_LOG("ERROR in 'LinearOperatorInverse::apply': Wrong storage format of Vectors. Aborting.\n");
				return false;
			}

			// copy d_nl as d
			domain_function_type& d = d_nl;

			UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN start rhs: \n" << d << " ----- END start rhs \n");
			UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN start sol: \n" << c_nl << " ----- END start sol \n");

			// build defect:  d := d_nl - J(u)*c_nl
			if(m_A->apply_sub(c_nl, d) != true)
			{
				UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to build defect. Aborting.\n");
				return false;
			}

			UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Start Defect: \n" << d << " ----- END Start Defect \n");

			// create correction, that has same memory pattern as u
			codomain_function_type c;
			if(c.clone_pattern(c_nl) != true)
			{
				UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to clone pattern for correction. Aborting.\n");
				return false;
			}

			// compute start norm ||d||_2
			number norm, norm_old, norm_start;
			norm = norm_old = norm_start = d.two_norm();

			// Print Start information
			if(m_verboseLevel >= 1) UG_LOG("\n    %%%%%%%%%% Iterative Linear Solver %%%%%%%%%%\n");
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

				// Compute a correction c := B*c using one the iterative step
				// Internally the defect is updated d := d - A*c = d - A*(x+c)
				if(m_iter->apply(d, c) == false)
				{
					UG_LOG("ERROR in 'LinearOperatorInverse::apply': Iterator Operator applied incorrectly. Aborting.\n");
					return false;
				}

				UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Correction of step " << i << ": \n" << c << " ----- END Correction of step " << i << "\n");
				UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Defect after step " << i << ": \n" << d << " ----- END Defect after step " << i << "\n");

				// add correction to solution
				c_nl += c;

				UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Sol after step " << i << ": \n" << c_nl << " ----- END Sol after step " << i << "\n");

				// compute global norm
				norm = d.two_norm();

				// print convergence rate
				if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// remember current norm
				norm_old = norm;
			}
			UG_ASSERT(0, "This line should never be reached.");
			return false;
		}

		// destructor
		virtual ~LinearSolver() {};

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

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__ */

/*
 * operator.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__OPERATOR__
#define __H__LIBDISCRETIZATION__OPERATOR__OPERATOR__

#include <iomanip>

#include "common/common.h"
#include "lib_discretization/assemble.h"

namespace ug{

// describes a mapping X->Y
template <typename X, typename Y>
class IDiscreteOperator
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(domain_function_type& u, codomain_function_type& d) = 0;

		// apply Operator, i.e. f := L(u);
		virtual bool apply(domain_function_type& u, codomain_function_type& d) = 0;

		// destructor
		virtual ~IDiscreteOperator() {};
};

// describes a mapping X->Y
template <typename X, typename Y>
class IDiscreteOperatorInverse
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init: This operator inverts the DiscreteOperator N: Y -> X
		virtual bool init(IDiscreteOperator<Y,X>& N) = 0;

		// prepare Operator
		virtual bool prepare(codomain_function_type& u) = 0;

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(codomain_function_type& u) = 0;

		// destructor
		virtual ~IDiscreteOperatorInverse() {};
};

// describes a mapping X->Y
template <typename X, typename Y>
class IDiscreteLinearizedOperator
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare J(u) for application of J(u)*c = d
		virtual bool prepare(domain_function_type& u, domain_function_type& c, codomain_function_type& d) = 0;

		// apply Operator, i.e. d = J(u)*c
		virtual bool apply(domain_function_type& c, codomain_function_type& d) = 0;

		// apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f) = 0;

		// destructor
		virtual ~IDiscreteLinearizedOperator() {};
};

/* This Operator type behaves different on application. It not only computes c = B(u)*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class IDiscreteLinearizedIteratorOperator
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// prepare for Linearized Operator J(u)
		virtual bool init(IDiscreteLinearizedOperator<Y, X>& J) = 0;

		// prepare B(u) for application of B(u)*d = c
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c) = 0;

		// compute new correction c = B(u)*d
		//    AND
		// update defect: d := d - J(u)*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c) = 0;

		// destructor
		virtual ~IDiscreteLinearizedIteratorOperator() {};
};


template <typename X, typename Y>
class IDiscreteLinearOperator : public IDiscreteLinearizedOperator<X,Y>
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// init
		virtual bool init() = 0;

		// prepare Operator
		virtual bool prepare(domain_function_type& u, codomain_function_type& f) = 0;

		// implement the interface for Linearized Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& c, codomain_function_type& d) {return prepare(c,d);};

		// apply Operator, i.e. f = L*u; (or d = J*c in case of Linearized Operator, i.e. u = c, f = d)
		virtual bool apply(domain_function_type& u, codomain_function_type& f) = 0;

		// apply Operator, i.e. f = f - L*u;
		virtual bool apply_sub(domain_function_type& u, codomain_function_type& f) = 0;

		// destructor
		virtual ~IDiscreteLinearOperator() {};
};

/* This Operator type behaves different on application. It not only computes c = B*d, but also changes d. */
/* It is used in iterative schemes. */
template <typename X, typename Y>
class IDiscreteLinearIteratorOperator : public IDiscreteLinearizedIteratorOperator<X,Y>
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		// prepare for Operator
		virtual bool init(IDiscreteLinearizedOperator<Y,X>& A) = 0;

		// prepare for correction and defect
		virtual bool prepare(domain_function_type& d, codomain_function_type& c) = 0;

		// Implement Interface for Linearized Operator
		virtual bool prepare(domain_function_type& u, domain_function_type& d, codomain_function_type& c){return prepare(d,c);}

		// compute new correction c = B*d
		//    AND
		// update defect: d := d - A*c
		virtual bool apply(domain_function_type& d, codomain_function_type& c) = 0;

		// destructor
		virtual ~IDiscreteLinearIteratorOperator() {};
};

template <typename X, typename Y>
class IDiscreteLinearizedOperatorInverse
{
	public:
		// export types:

		// domain space
		typedef X domain_function_type;

		// range space
		typedef Y codomain_function_type;

	public:
		IDiscreteLinearizedOperatorInverse( IDiscreteLinearizedIteratorOperator<X,Y>& B, int maxIter, number absTol, number relTol,
										int verboseLevel = 0) :
			m_verboseLevel(verboseLevel), m_iter(&B),
			m_maxIter(maxIter), m_absTol(absTol), m_relTol(relTol)
			{};

		virtual bool init(IDiscreteLinearizedOperator<Y,X>& A)
		{
			m_A = &A;
			return true;
		}

		// prepare Operator
		virtual bool prepare(codomain_function_type& u, domain_function_type& d, codomain_function_type& c)
		{
			// TODO: Do we assume, that operator has been prepared? Do we have to prepare it here?

			// init iterator B for operator A
			if(m_iter->init(*m_A) != true)
			{
				UG_LOG("ERROR in 'DiscreteLinearizedOperatorInverse::prepare': Cannot init Iterator Operator for Operator A.\n");
				return false;
			}

			// prepare iterator B for d and c
			if(m_iter->prepare(u, d, c) != true)
			{
				UG_LOG("ERROR in 'DiscreteLinearizedOperatorInverse::prepare': Cannot prepare Iterator Operator.\n");
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
			// copy d_nl as d
			domain_function_type& d = d_nl;

			UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN start rhs: \n" << d << " ----- END start rhs \n");
			UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN start sol: \n" << c_nl << " ----- END start sol \n");

			// build defect:  d := d_nl - J(u)*c_nl
			m_A->apply_sub(c_nl, d);

			UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Start Defect: \n" << d << " ----- END Start Defect \n");

			// create correction, that has same memory pattern as u
			codomain_function_type c(c_nl);

#ifdef UG_PARALLEL
			// make defect unique
			d.parallel_additive_to_unique();
#endif

			// compute start norm ||d||_2
			number norm, norm_old, norm_start;
			norm = norm_old = norm_start = d.norm();

			// Print Start information
			if(m_verboseLevel >= 1) UG_LOG("\n    %%%%%%%%%% Iterative Linear Solver %%%%%%%%%%\n");
			if(m_verboseLevel >= 2) UG_LOG("    %   Iter     Defect         Rate \n");
			if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << 0 << ":  " << std::scientific << norm_old <<  "      -------" << std::endl);

			// Iteration loop
			for(int i = 1; ; ++i)
			{
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

				// check that defect is a still a valid number
				if(!is_valid_number(norm))
				{
					if(m_verboseLevel >= 1) UG_LOG("    %%%%% Defect " << norm << " is not a valid number. Linear Solver did NOT CONVERGE. %%%%%\n\n");
					return false;
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
					UG_LOG("ERROR in 'DiscreteLinearOperatorInverse::apply': Iterator Operator applied incorrectly. Aborting.\n");
					return false;
				}

				UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Correction of step " << i << ": \n" << c << " ----- END Correction of step " << i << "\n");
				UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Defect after step " << i << ": \n" << d << " ----- END Defect after step " << i << "\n");

				// add correction to solution
				c_nl += c;

				UG_DLOG(LIB_DISC_OPERATOR_INVERSE, 10, " ----- BEGIN Sol after step " << i << ": \n" << c_nl << " ----- END Sol after step " << i << "\n");

#ifdef UG_PARALLEL
				// make defect unique
				d.parallel_additive_to_unique();
#endif

				// compute new defect norm
				double tNormLocal = (double)d.norm();
				double tNormGlobal;
				pcl::AllReduce(&tNormLocal, &tNormGlobal, 1, PCL_DT_DOUBLE, PCL_RO_SUM);
				norm = (number)tNormGlobal;

				// print convergence rate
				if(m_verboseLevel >= 2) UG_LOG("    % " << std::setw(4) << i << ":  " << std::scientific << norm << "    " << norm/norm_old << std::endl);

				// remember current norm
				norm_old = norm;
			}
			UG_ASSERT(0, "This line should never be reached.");
			return false;
		}

		// destructor
		virtual ~IDiscreteLinearizedOperatorInverse() {};

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

			return value >= std::numeric_limits<number>::min() && value <= std::numeric_limits<number>::max() && value == value;
		}

		// Discribes, how many output is printed. (0 = nothing, 1 = major informations, 2 = all)
		int m_verboseLevel;

	protected:
		// Operator that is inverted by this Inverse Operator
		IDiscreteLinearizedOperator<Y,X>* m_A;

		// Iterator used in the iterative scheme to compute the correction and update the defect
		IDiscreteLinearizedIteratorOperator<X,Y>* m_iter;

		// maximal number of iterations
		int m_maxIter;

		// absolute defect to be reached
		number m_absTol;

		// relative defect to be reached
		number m_relTol;

};

}

#include "linear_operator/interpolation_operator.h"
#include "linear_operator/projection_operator.h"
#include "linear_operator/transfer_operator.h"
#include "linear_operator/assembled_linear_operator.h"
#include "linear_operator/multi_grid_solver/mg_solver.h"

#include "non_linear_operator/assembled_non_linear_operator.h"
#include "non_linear_operator/newton_solver/newton.h"

#endif /* __H__LIBDISCRETIZATION__OPERATOR__OPERATOR__ */

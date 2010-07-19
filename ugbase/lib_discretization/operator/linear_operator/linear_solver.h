/*
 * linear_solver.h
 *
 *  Created on: 22.02.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__
#define __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__

#include "lib_discretization/operator/operator.h"

namespace ug{

template <typename TFunction>
class LinearSolver : public ILinearizedOperatorInverse<TFunction, TFunction>
{
	public:
		// domain space
		typedef TFunction domain_function_type;

		// range space
		typedef TFunction codomain_function_type;

	public:
		LinearSolver( 	ILinearizedIteratorOperator<TFunction,TFunction>* Precond,
						ConvergenceCheck<TFunction>& ConvCheck) :
							m_pPrecond(Precond), m_ConvCheck(ConvCheck)
			{};

		virtual bool init(ILinearizedOperator<TFunction, TFunction>& A)
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
			if(m_pPrecond != NULL)
			{
				if(m_pPrecond->init(*m_A) != true)
				{UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot init "
							"Iterator Operator for Operator A.\n");return false;}

				// prepare iterator B for d_nl and c_nl
				if(m_pPrecond->prepare(u, d_nl, c_nl) != true)
				{UG_LOG("ERROR in 'LinearizedOperatorInverse::prepare': Cannot "
							"prepare Iterator Operator.\n"); return false;}
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
			if(!d_nl.has_storage_type(PST_ADDITIVE))
				{UG_LOG("ERROR in 'LinearSolver::apply': d must be additive. Aborting.\n"); return false;}
			if(!c_nl.has_storage_type(PST_CONSISTENT))
				{UG_LOG("ERROR in 'LinearSolver::apply': c must be consistent. Aborting.\n"); return false;}

			// copy d_nl as d
			domain_function_type& d = d_nl;

			// build defect:  d := d_nl - J(u)*c_nl
			if(!m_A->apply_sub(c_nl, d))
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to build defect. Aborting.\n"); return false;}

			// create correction, that has same memory pattern as u
			codomain_function_type c;
			if(!c.clone_pattern(c_nl))
				{UG_LOG("ERROR in 'LinearOperatorInverse::apply': Unable to "
							"clone pattern for correction. Aborting.\n"); return false;}

			m_ConvCheck.set_offset(3);
			m_ConvCheck.set_symbol('%');
			m_ConvCheck.set_name("Iterative Linear Solver");
			m_ConvCheck.start(d);

			// Iteration loop
			while(!m_ConvCheck.iteration_ended())
			{
				// Compute a correction c := B*c using one the iterative step
				// Internally the defect is updated d := d - A*c = d - A*(x+c)
				if(m_pPrecond != NULL)
					if(!m_pPrecond->apply(d, c, true))
						{UG_LOG("ERROR in 'LinearOperatorInverse::apply': Iterator Operator "
									"applied incorrectly. Aborting.\n"); return false;}

				// add correction to solution
				c_nl += c;

				// check convergence
				m_ConvCheck.update(d);
			}

			return m_ConvCheck.post();
		}

		// destructor
		virtual ~LinearSolver() {};

	protected:
		// Operator that is inverted by this Inverse Operator
		ILinearizedOperator<TFunction,TFunction>* m_A;

		// Iterator used in the iterative scheme to compute the correction and update the defect
		ILinearizedIteratorOperator<TFunction,TFunction>* m_pPrecond;

		// Convergence Check
		ConvergenceCheck<TFunction>& m_ConvCheck;
};

} // end namespace ug

#endif /* __H__LIBDISCRETIZATION__OPERATOR__LINEAR_OPERATOR__LINEAR_SOLVER__ */

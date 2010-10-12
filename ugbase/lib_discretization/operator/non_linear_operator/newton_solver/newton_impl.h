/*
 * newton_impl.h
 *
 *  Created on: 18.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__
#define __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__

#include "newton.h"
#include "lib_discretization/io/vtkoutput.h"
#include "lib_discretization/function_spaces/grid_function_util.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra>
bool
NewtonSolver<TDoFDistribution, TAlgebra>::
init(IOperator<vector_type, vector_type>& N)
{
	m_N = dynamic_cast<AssembledOperator<dof_distribution_type, algebra_type>* >(&N);
	if(m_N == NULL)
	{
		UG_LOG("NewtonSolver currently only works for AssembledDiscreteOperator.\n");
		return false;
	}

	m_pAss = m_N->get_assemble();
	return true;
}


template <typename TDoFDistribution, typename TAlgebra>
bool
NewtonSolver<TDoFDistribution, TAlgebra>::
allocate_memory(const vector_type& u)
{
	// Jacobian
	m_J = new AssembledLinearOperator<dof_distribution_type, algebra_type>(*m_pAss);
	m_J->set_dof_distribution(*m_N->get_dof_distribution());

	// defect
	m_d.resize(u.size()); m_d = u;

	// correction
	m_c.resize(u.size()); m_c = u;

	if(m_J == NULL)
	{
		UG_ASSERT(0, "Cannot allocate memory.");
		return false;
	}

	m_allocated = true;
	return true;
}

template <typename TDoFDistribution, typename TAlgebra>
bool
NewtonSolver<TDoFDistribution, TAlgebra>::
deallocate_memory()
{
	if(m_allocated)
	{
		m_d.destroy();
		m_c.destroy();
		delete m_J;
	}
	m_allocated = false;
	return true;
}


template <typename TDoFDistribution, typename TAlgebra>
bool
NewtonSolver<TDoFDistribution, TAlgebra>::
prepare(vector_type& u)
{
	if(!m_allocated)
	{
		if(allocate_memory(u) != true)
		{
			UG_LOG("NewtonSolver: Cannot allocate memory.\n");
			return false;
		}
	}

	return true;
}


template <typename TDoFDistribution, typename TAlgebra>
NewtonSolver<TDoFDistribution, TAlgebra>::
~NewtonSolver()
{
	if(m_allocated)
	{
		if(!deallocate_memory())
			UG_ASSERT(0, "Cannot deallocate memory");
	}
}


template <typename TDoFDistribution, typename TAlgebra>
bool
NewtonSolver<TDoFDistribution, TAlgebra>::
apply(vector_type& u)
{
	if(m_pLinearSolver == NULL)
	{
		UG_LOG("ERROR in 'NewtonSolver::apply': Linear Solver not set.\n");
		return false;
	}
	if(m_pConvCheck == NULL)
	{
		UG_LOG("ERROR in 'NewtonSolver::apply': Convergence Check not set.\n");
		return false;
	}

	// Compute first Defect
	if(m_N->prepare(m_d, u) != true)
		{UG_LOG("NewtonSolver::apply: Cannot prepare Non-linear Operator.\n"); return false;}
	if(m_N->apply(m_d, u) != true)
		{UG_LOG("NewtonSolver::apply: Cannot apply Non-linear Operator to compute start defect.\n"); return false;}

	// start convergence check
	m_pConvCheck->set_offset(3);
	m_pConvCheck->set_symbol('#');
	m_pConvCheck->set_name("Newton Solver");
	m_pConvCheck->start(m_d);

	// copy pattern
	vector_type s; s.resize(u.size()); s = u;

	//loop iteration
	while(!m_pConvCheck->iteration_ended())
	{
		// set c = 0
		if(!m_c.set(0.0))
			{UG_LOG("NewtonSolver::apply: Cannot reset correction to zero.\n"); return false;}

		// Compute Jacobian
		if(!m_J->init(u))
			{UG_LOG("NewtonSolver::apply: Cannot prepare Jacobi Operator.\n"); return false;}

		// Init Jacobi Inverse
		if(!m_pLinearSolver->init(*m_J, u))
			{UG_LOG("NewtonSolver::apply: Cannot init Inverse Linear "
					"Operator for Jacobi-Operator.\n"); return false;}

		// Solve Linearized System
		if(!m_pLinearSolver->apply(m_c, m_d))
			{UG_LOG("NewtonSolver::apply: Cannot apply Inverse Linear "
					"Operator for Jacobi-Operator.\n"); return false;}

		// Line Search
		if(m_pLineSearch != NULL)
		{
			m_pLineSearch->set_offset("   #  ");
			if(!m_pLineSearch->search(*m_N, u, m_c, m_d, m_pConvCheck->defect()))
				{UG_LOG("Newton Solver did not converge.\n"); return false;}
		}
		// no line search: Compute new defect
		else
		{
			// update solution
			u -= m_c;

			// compute new Defect
			if(!m_N->prepare(m_d, u))
				{UG_LOG("NewtonSolver::apply: Cannot prepare Non-linear"
						" Operator for defect computation.\n"); return false;}
			if(!m_N->apply(m_d, u))
				{UG_LOG("NewtonSolver::apply: Cannot apply Non-linear Operator "
						"to compute defect.\n"); return false;}
		}

		// check convergence
		m_pConvCheck->update(m_d);
	}

	return m_pConvCheck->post();
}

}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__ */


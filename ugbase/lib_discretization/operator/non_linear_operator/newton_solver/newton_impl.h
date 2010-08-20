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

template <typename TFunction>
bool
NewtonSolver<TFunction>::
init(IOperator<function_type, function_type>& N)
{
	m_N = dynamic_cast<AssembledOperator<function_type>* >(&N);
	if(m_N == NULL)
	{
		UG_LOG("NewtonSolver currently only works for AssembledDiscreteOperator.\n");
		return false;
	}
	m_ass = m_N->get_assemble();
	return true;
}


template <typename TFunction>
bool
NewtonSolver<TFunction>::
allocate_memory(const function_type& u)
{
	// Jacobian
	m_J = new AssembledLinearizedOperator<function_type>(*m_ass);

	// defect
	m_d = new function_type;
	m_d->clone_pattern(u);

	// correction
	m_c = new function_type;
	m_c->clone_pattern(u);

	if(m_d == NULL || m_c == NULL || m_J == NULL)
	{
		UG_ASSERT(0, "Cannot allocate memory.");
		return false;
	}

	m_allocated = true;
	return true;
}

template <typename TFunction>
bool
NewtonSolver<TFunction>::
deallocate_memory()
{
	if(m_allocated)
	{
		delete m_d; delete m_c;
		delete m_J;
	}
	m_allocated = false;
	return true;
}


template <typename TFunction>
bool
NewtonSolver<TFunction>::
prepare(function_type& u)
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


template <typename TFunction>
NewtonSolver<TFunction>::
~NewtonSolver()
{
	if(m_allocated)
	{
		if(!deallocate_memory())
			UG_ASSERT(0, "Cannot deallocate memory");
	}
}


template <typename TFunction>
bool
NewtonSolver<TFunction>::
apply(function_type& u)
{
	// Compute first Defect
	if(m_N->prepare(u, *m_d) != true)
		{UG_LOG("NewtonSolver::apply: Cannot prepare Non-linear Operator.\n"); return false;}
	if(m_N->apply(u, *m_d) != true)
		{UG_LOG("NewtonSolver::apply: Cannot apply Non-linear Operator to compute start defect.\n"); return false;}

	// start convergence check
	m_ConvCheck.set_offset(3);
	m_ConvCheck.set_symbol('#');
	m_ConvCheck.set_name("Newton Solver");
	m_ConvCheck.start(*m_d);

	function_type s;
	s.clone_pattern(u);

	//loop iteration
	while(!m_ConvCheck.iteration_ended())
	{
		// set c = 0
		if(!m_c->set(0.0))
			{UG_LOG("NewtonSolver::apply: Cannot reset correction to zero.\n"); return false;}

		// Compute Jacobian
		if(!m_J->prepare(u, *m_c, *m_d))
			{UG_LOG("NewtonSolver::apply: Cannot prepare Jacobi Operator.\n"); return false;}

		WriteMatrixToConnectionViewer("NewtonMat.mat", m_J->get_matrix(), u);

		// Init Jacobi Inverse
		if(!m_LinearSolver.init(*m_J))
			{UG_LOG("NewtonSolver::apply: Cannot init Inverse Linear "
					"Operator for Jacobi-Operator.\n"); return false;}
		if(!m_LinearSolver.prepare(u, *m_d, *m_c))
			{UG_LOG("NewtonSolver::apply: Cannot prepare Inverse Linear "
					"Operator for Jacobi-Operator.\n"); return false;}

		// Solve Linearized System
		if(!m_LinearSolver.apply(*m_d, *m_c))
			{UG_LOG("NewtonSolver::apply: Cannot apply Inverse Linear "
					"Operator for Jacobi-Operator.\n"); return false;}

		// Line Search
		if(m_LineSearch != NULL)
		{
			m_LineSearch->set_offset("   #  ");
			if(!m_LineSearch->search(*m_N, u, *m_c, *m_d, m_ConvCheck.defect()))
				{UG_LOG("Newton Solver did not converge.\n"); return false;}
		}
		// no line search: Compute new defect
		else
		{
			// update solution
			u -= *m_c;

			// compute new Defect
			if(!m_N->prepare(u, *m_d))
				{UG_LOG("NewtonSolver::apply: Cannot prepare Non-linear"
						" Operator for defect computation.\n"); return false;}
			if(!m_N->apply(u, *m_d))
				{UG_LOG("NewtonSolver::apply: Cannot apply Non-linear Operator "
						"to compute defect.\n"); return false;}
		}

		// check convergence
		m_ConvCheck.update(*m_d);
	}

	return m_ConvCheck.post();
}

}

#endif /* __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__ */


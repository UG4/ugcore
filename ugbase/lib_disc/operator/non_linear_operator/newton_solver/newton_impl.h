/*
 * newton_impl.h
 *
 *  Created on: 18.06.2010
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__

#include <iostream>
#include <sstream>

#include "newton.h"
#include "lib_disc/function_spaces/grid_function_util.h"

#define PROFILE_NEWTON
#ifdef PROFILE_NEWTON
	#define NEWTON_PROFILE_FUNC()		PROFILE_FUNC()
	#define NEWTON_PROFILE_BEGIN(name)	PROFILE_BEGIN(name)
	#define NEWTON_PROFILE_END()		PROFILE_END()
#else
	#define NEWTON_PROFILE_FUNC()
	#define NEWTON_PROFILE_BEGIN(name)
	#define NEWTON_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
bool
NewtonSolver<TAlgebra>::
init(IOperator<vector_type, vector_type>& N)
{
	m_N = dynamic_cast<AssembledOperator<TAlgebra>* >(&N);
	if(m_N == NULL)
	{
		UG_LOG("NewtonSolver currently only works for AssembledDiscreteOperator.\n");
		return false;
	}

	m_pAss = m_N->get_assemble();
	return true;
}


template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::allocate_memory(const vector_type& u)
{
	// Jacobian
	m_J = new AssembledLinearOperator<TAlgebra>(*m_pAss);
	m_J->set_level(m_N->level());

	// defect
	m_d.resize(u.size()); m_d = u;

	// correction
	m_c.resize(u.size()); m_c = u;

	if(m_J == NULL)
		UG_THROW_FATAL("Cannot allocate memory.");

	m_allocated = true;
	return true;
}

template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::deallocate_memory()
{
	if(m_allocated)
		delete m_J;
	m_allocated = false;
	return true;
}


template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::prepare(vector_type& u)
{
	if(!m_allocated)
	{
		if(allocate_memory(u) != true)
		{
			UG_LOG("NewtonSolver: Cannot allocate memory.\n");
			return false;
		}
	}

//	Check for linear solver
	if(m_pLinearSolver == NULL)
	{
		UG_LOG("ERROR in 'NewtonSolver::prepare': Linear Solver not set.\n");
		return false;
	}

//	Check if ConvCheck has been set
	if(m_pConvCheck == NULL)
	{
		UG_LOG("ERROR in 'NewtonSolver::prepare': Convergence Check not set.\n");
		return false;
	}

//	Set dirichlet values
	m_N->prepare(m_d, u);

	return true;
}


template <typename TAlgebra>
NewtonSolver<TAlgebra>::~NewtonSolver()
{
	if(m_allocated)
	{
		if(!deallocate_memory())
			UG_ASSERT(0, "Cannot deallocate memory");
	}

}


template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::apply(vector_type& u)
{
//	increase call count
	m_dgbCall++;

//	resize
	m_d.resize(u.size());
	m_c.resize(u.size());

// 	Compute first Defect
	NEWTON_PROFILE_BEGIN(NewtonComputeDefect1);
	m_N->apply(m_d, u);
	NEWTON_PROFILE_END();

//	write start defect for debug
	int loopCnt = 0;
	char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
	std::string name("NEWTON_Defect");
	name.append(ext);
	write_debug(m_d, name.c_str());
	write_debug(u, "NEWTON_StartSolution");

// 	increase offset of output for linear solver
	IConvergenceCheck* pLinConvCheck = m_pLinearSolver->get_convergence_check();
	int iLinSolverOffset = 0;
	if(pLinConvCheck != NULL)
	{
		iLinSolverOffset = pLinConvCheck->get_offset();
		pLinConvCheck->set_offset( m_pConvCheck->get_offset() + 3);
	}

// 	set info string indicating the used linear solver
	std::stringstream ss; ss << "(Linear Solver: " << m_pLinearSolver->name() << ")";
	m_pConvCheck->set_info(ss.str());

// 	copy pattern
	vector_type s; s.resize(u.size()); s = u;

// 	start convergence check
	m_pConvCheck->start(m_d);

//	loop iteration
	while(!m_pConvCheck->iteration_ended())
	{
	// 	set c = 0
		NEWTON_PROFILE_BEGIN(NewtonSetCorretionZero);
		if(!m_c.set(0.0))
		{
			UG_LOG("ERROR in 'NewtonSolver::apply':"
					" Cannot reset correction to zero.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	// 	Compute Jacobian
		NEWTON_PROFILE_BEGIN(NewtonComputeJacobian);
		m_J->init(u);
		NEWTON_PROFILE_END();

	//	Write Jacobian for debug
		std::string matname("NEWTON_Jacobian");
		matname.append(ext);
		write_debug(m_J->get_matrix(), matname.c_str());

	// 	Init Jacobi Inverse
		NEWTON_PROFILE_BEGIN(NewtonPrepareLinSolver);
		if(!m_pLinearSolver->init(*m_J, u))
		{
			UG_LOG("ERROR in 'NewtonSolver::apply': Cannot init Inverse Linear "
					"Operator for Jacobi-Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	// 	Solve Linearized System
		NEWTON_PROFILE_BEGIN(NewtonApplyLinSolver);
		if(!m_pLinearSolver->apply(m_c, m_d))
		{
			UG_LOG("ERROR in 'NewtonSolver::apply': Cannot apply Inverse Linear "
					"Operator for Jacobi-Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();

	//	store convergence history
		IConvergenceCheck* pLinConvCheck = m_pLinearSolver->get_convergence_check();
		const int numSteps = pLinConvCheck->step();
		if(loopCnt >= (int)m_vTotalLinSolverSteps.size()) m_vTotalLinSolverSteps.resize(loopCnt+1);
		if(loopCnt >= (int)m_vLinSolverCalls.size()) m_vLinSolverCalls.resize(loopCnt+1, 0);
		m_vTotalLinSolverSteps[loopCnt] += numSteps;
		m_vLinSolverCalls[loopCnt] += 1;

	// 	Line Search
		if(m_pLineSearch != NULL)
		{
			m_pLineSearch->set_offset("   #  ");
			NEWTON_PROFILE_BEGIN(NewtonLineSearch);
			if(!m_pLineSearch->search(*m_N, u, m_c, m_d, m_pConvCheck->defect()))
			{
				UG_LOG("ERROR in 'NewtonSolver::apply': "
						"Newton Solver did not converge.\n");
				return false;
			}
			NEWTON_PROFILE_END();
		}
	// 	No line search: Compute new defect
		else
		{
		// 	update solution
			u -= m_c;

		// 	compute new Defect
			NEWTON_PROFILE_BEGIN(NewtonComputeDefect);
			m_N->prepare(m_d, u);
			m_N->apply(m_d, u);
			NEWTON_PROFILE_END();
		}

	//	update counter
		loopCnt++;
		sprintf(ext, "_iter%03d", loopCnt);

	//	write defect for debug
		std::string name("NEWTON_Defect"); name.append(ext);
		write_debug(m_d, name.c_str());
		std::string name2("NEWTON_Correction"); name2.append(ext);
		write_debug(m_c, name2.c_str());

	// 	check convergence
		m_pConvCheck->update(m_d);
	}

	// reset offset of output for linear solver to previous value
	if(pLinConvCheck != NULL)
	{
		pLinConvCheck->set_offset(iLinSolverOffset);
	}

	return m_pConvCheck->post();
}

template <typename TAlgebra>
void
NewtonSolver<TAlgebra>::
print_average_convergence() const
{
	UG_LOG("\nNewton solver convergence history:\n");
	UG_LOG("Newton Step | Num Calls | Total Lin Iters | Avg Lin Iters \n");
	int allCalls = 0, allSteps = 0;
	for(int call = 0; call < (int)m_vLinSolverCalls.size(); ++call)
	{
		UG_LOG( " " << std::setw(10) << call+1 << " | ");
		UG_LOG(std::setw(9) << m_vLinSolverCalls[call] << " | ");
		allCalls += m_vLinSolverCalls[call];
		UG_LOG(std::setw(15) << m_vTotalLinSolverSteps[call] << " | ");
		allSteps += m_vTotalLinSolverSteps[call];
		UG_LOG(std::setw(13) << std::setprecision(2) << std::fixed << m_vTotalLinSolverSteps[call] / (double)m_vLinSolverCalls[call]);
		UG_LOG("\n");
	}
	UG_LOG( "        all | ");
	UG_LOG(std::setw(9) << allCalls << " | ");
	UG_LOG(std::setw(15) << allSteps << " | ");
	UG_LOG(std::setw(13) << std::setprecision(2) << std::fixed << allSteps / (double)allCalls);
	UG_LOG("\n");
}

template <typename TAlgebra>
void
NewtonSolver<TAlgebra>::
clear_average_convergence()
{
	m_vLinSolverCalls.clear();
	m_vTotalLinSolverSteps.clear();
}


}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__ */


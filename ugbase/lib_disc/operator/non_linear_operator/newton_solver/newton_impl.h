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
#include "common/util/string_util.h"

#define PROFILE_NEWTON
#ifdef PROFILE_NEWTON
	#define NEWTON_PROFILE_FUNC()		PROFILE_FUNC_GROUP("Newton")
	#define NEWTON_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "Newton")
	#define NEWTON_PROFILE_END()		PROFILE_END()
#else
	#define NEWTON_PROFILE_FUNC()
	#define NEWTON_PROFILE_BEGIN(name)
	#define NEWTON_PROFILE_END()
#endif

namespace ug{

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver,
			SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
			SmartPtr<ILineSearch<vector_type> > spLineSearch) :
			m_spLinearSolver(LinearSolver),
			m_spConvCheck(spConvCheck),
			m_spLineSearch(spLineSearch),
			m_N(NULL),
			m_J(NULL),
			m_spAss(NULL),
			m_dgbCall(0)
{};

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver() :
	m_spLinearSolver(NULL),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_spLineSearch(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(NULL),
	m_dgbCall(0)
{};

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver(SmartPtr<IOperator<vector_type> > N) :
	m_spLinearSolver(NULL),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_spLineSearch(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(NULL),
	m_dgbCall(0)
{
	init(N);
};

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver(SmartPtr<IAssemble<TAlgebra> > spAss) :
	m_spLinearSolver(NULL),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_spLineSearch(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(NULL),
	m_dgbCall(0)
{
	m_spAss = spAss;
	m_N = SmartPtr<AssembledOperator<TAlgebra> >(new AssembledOperator<TAlgebra>(m_spAss));
};

template <typename TAlgebra>
void NewtonSolver<TAlgebra>::
set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck)
{
	m_spConvCheck = spConvCheck;
	m_spConvCheck->set_offset(3);
	m_spConvCheck->set_symbol('#');
	m_spConvCheck->set_name("Newton Solver");
}

template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::init(SmartPtr<IOperator<vector_type> > N)
{
	NEWTON_PROFILE_BEGIN(NewtonSolver_init);
	m_N = N.template cast_dynamic<AssembledOperator<TAlgebra> >();
	if(m_N.invalid())
		UG_THROW("NewtonSolver: currently only works for AssembledDiscreteOperator.");

	m_spAss = m_N->discretization();
	return true;
}

template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::prepare(vector_type& u)
{
//	todo: maybe remove this from interface
	return true;
}

template <typename TAlgebra>
bool NewtonSolver<TAlgebra>::apply(vector_type& u)
{
	NEWTON_PROFILE_BEGIN(NewtonSolver_apply);
//	increase call count
	m_dgbCall++;

//	Check for linear solver
	if(m_spLinearSolver.invalid())
		UG_THROW("NewtonSolver::apply: Linear Solver not set.");

//	Jacobian
	if(m_J.invalid() || m_J->discretization() != m_spAss) {
		m_J = make_sp(new AssembledLinearOperator<TAlgebra>(m_spAss));
	}
	m_J->set_level(m_N->level());

//	create tmp vectors
	SmartPtr<vector_type> spD = u.clone_without_values();
	SmartPtr<vector_type> spC = u.clone_without_values();

//	Set dirichlet values
	try{
		m_N->prepare(u);
	}
	UG_CATCH_THROW("NewtonSolver::prepare: Prepare of Operator failed.");

// 	Compute first Defect
	try{
	NEWTON_PROFILE_BEGIN(NewtonComputeDefect1);
	m_N->apply(*spD, u);
	NEWTON_PROFILE_END();
	}UG_CATCH_THROW("NewtonSolver::apply: Computation of Start-Defect failed.");

//	write start defect for debug
	int loopCnt = 0;
	char ext[20]; sprintf(ext, "_iter%03d", loopCnt);
	std::string name("NEWTON_Defect");
	name.append(ext);
	write_debug(*spD, name.c_str());
	write_debug(u, "NEWTON_StartSolution");

// 	increase offset of output for linear solver
	const int stdLinOffset = m_spLinearSolver->standard_offset();
	m_spLinearSolver->convergence_check()->set_offset(stdLinOffset + 3);

// 	set info string indicating the used linear solver
	std::stringstream ss; ss << "(Linear Solver: " << m_spLinearSolver->name() << ")";
	m_spConvCheck->set_info(ss.str());

// 	start convergence check
	m_spConvCheck->start(*spD);

	for(size_t i = 0; i < m_stepUpdate.size(); ++i)
		m_stepUpdate[i]->update();

//	loop iteration
	while(!m_spConvCheck->iteration_ended())
	{
	// 	set c = 0
		NEWTON_PROFILE_BEGIN(NewtonSetCorretionZero);
		spC->set(0.0);
		NEWTON_PROFILE_END();

		for(size_t i = 0; i < m_innerStepUpdate.size(); ++i)
			m_innerStepUpdate[i]->update();

	// 	Compute Jacobian
		try{
		NEWTON_PROFILE_BEGIN(NewtonComputeJacobian);
		m_J->init(u);
		NEWTON_PROFILE_END();
		}UG_CATCH_THROW("NewtonSolver::apply: Initialization of Jacobian failed.");

	//	Write Jacobian for debug
		std::string matname("NEWTON_Jacobian");
		matname.append(ext);
		write_debug(m_J->get_matrix(), matname.c_str());

	// 	Init Jacobi Inverse
		try{
		NEWTON_PROFILE_BEGIN(NewtonPrepareLinSolver);
		if(!m_spLinearSolver->init(m_J, u))
		{
			UG_LOG("ERROR in 'NewtonSolver::apply': Cannot init Inverse Linear "
					"Operator for Jacobi-Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();
		}UG_CATCH_THROW("NewtonSolver::apply: Initialization of Linear Solver failed.");

	// 	Solve Linearized System
		try{
		NEWTON_PROFILE_BEGIN(NewtonApplyLinSolver);
		if(!m_spLinearSolver->apply(*spC, *spD))
		{
			UG_LOG("ERROR in 'NewtonSolver::apply': Cannot apply Inverse Linear "
					"Operator for Jacobi-Operator.\n");
			return false;
		}
		NEWTON_PROFILE_END();
		}UG_CATCH_THROW("NewtonSolver::apply: Application of Linear Solver failed.");

	//	store convergence history
		const int numSteps = m_spLinearSolver->step();
		if(loopCnt >= (int)m_vTotalLinSolverSteps.size()) m_vTotalLinSolverSteps.resize(loopCnt+1);
		if(loopCnt >= (int)m_vLinSolverCalls.size()) m_vLinSolverCalls.resize(loopCnt+1, 0);
		if(loopCnt >= (int)m_vLinSolverRates.size()) m_vLinSolverRates.resize(loopCnt+1, 0);
		m_vTotalLinSolverSteps[loopCnt] += numSteps;
		m_vLinSolverCalls[loopCnt] += 1;
		m_vLinSolverRates[loopCnt] += m_spLinearSolver->convergence_check()->avg_rate();

	// 	Line Search
		try{
		if(m_spLineSearch.valid())
		{
			m_spLineSearch->set_offset("   #  ");
			NEWTON_PROFILE_BEGIN(NewtonLineSearch);
			if(!m_spLineSearch->search(m_N, u, *spC, *spD, m_spConvCheck->defect()))
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
			u -= *spC;

		// 	compute new Defect
			NEWTON_PROFILE_BEGIN(NewtonComputeDefect);
			m_N->prepare(u);
			m_N->apply(*spD, u);
			NEWTON_PROFILE_END();
		}
		}UG_CATCH_THROW("NewtonSolver::apply: Line Search update failed.");


	//	update counter
		loopCnt++;
		sprintf(ext, "_iter%03d", loopCnt);

	// 	check convergence
		m_spConvCheck->update(*spD);
		if(loopCnt-1 >= (int)m_vNonLinSolverRates.size()) m_vNonLinSolverRates.resize(loopCnt, 0);
		m_vNonLinSolverRates[loopCnt-1] += m_spConvCheck->rate();

	//	write defect for debug
		std::string name("NEWTON_Defect"); name.append(ext);
		write_debug(*spD, name.c_str());
		std::string name2("NEWTON_Correction"); name2.append(ext);
		write_debug(*spC, name2.c_str());
		std::string name3("NEWTON_Solution"); name3.append(ext);
		write_debug(u, name3.c_str());
	}

	// reset offset of output for linear solver to previous value
	m_spLinearSolver->convergence_check()->set_offset(stdLinOffset);

	return m_spConvCheck->post();
}

template <typename TAlgebra>
void NewtonSolver<TAlgebra>::print_average_convergence() const
{
	UG_LOG("\nNewton solver convergence history:\n");
	UG_LOG("Newton Step | Num Calls | Total Lin Iters | Avg Lin Iters | Avg Nonlin Rates | Avg Lin Rates \n");
	int allCalls = 0, allLinSteps = 0;
	number allLinRatesProduct = 1.0, allNonLinRatesProduct = 1.0;
	for(int call = 0; call < (int)m_vLinSolverCalls.size(); ++call)
	{
		UG_LOG( " " << std::setw(10) << call+1 << " | ");
		UG_LOG(std::setw(9) << m_vLinSolverCalls[call] << " | ");
		allCalls += m_vLinSolverCalls[call];
		UG_LOG(std::setw(15) << m_vTotalLinSolverSteps[call] << " | ");
		allLinSteps += m_vTotalLinSolverSteps[call];
		UG_LOG(std::setw(13) << std::setprecision(2) << std::fixed << m_vTotalLinSolverSteps[call] / (double)m_vLinSolverCalls[call] << " | ");
		allNonLinRatesProduct *= pow((number)m_vNonLinSolverRates[call]/(double)m_vLinSolverCalls[call],(double)m_vLinSolverCalls[call]);
		UG_LOG(std::setw(16) << std::setprecision(6) << std::scientific << m_vNonLinSolverRates[call] / (double)m_vLinSolverCalls[call] << " | ");
		allLinRatesProduct *= (number)std::pow((number)m_vLinSolverRates[call]/(double)m_vLinSolverCalls[call],(number)m_vTotalLinSolverSteps[call]);
		UG_LOG(std::setw(13) << std::setprecision(6) << std::scientific << m_vLinSolverRates[call] / (double)m_vLinSolverCalls[call]);
		UG_LOG("\n");
	}
	UG_LOG( "        all | ");
	UG_LOG(std::setw(9) << allCalls << " | ");
	UG_LOG(std::setw(15) << allLinSteps << " | ");
	UG_LOG(std::setw(13) << std::setprecision(2) << std::fixed << allLinSteps / (number)allCalls << " | ");
	UG_LOG(std::setw(16) << std::setprecision(6) << std::scientific << std::pow((number)allNonLinRatesProduct,(number)1.0/(number)allCalls) << " | ");
	UG_LOG(std::setw(13) << std::setprecision(6) << std::scientific << std::pow((number)allLinRatesProduct,(number)1.0/(number)allLinSteps));
	UG_LOG("\n");
}

template <typename TAlgebra>
size_t NewtonSolver<TAlgebra>::num_newton_steps() const
{
	return m_vLinSolverCalls.size();
}

template <typename TAlgebra>
int NewtonSolver<TAlgebra>::num_linsolver_calls(size_t call) const
{
	return m_vLinSolverCalls[call];
}

template <typename TAlgebra>
int NewtonSolver<TAlgebra>::num_linsolver_steps(size_t call) const
{
	return m_vTotalLinSolverSteps[call];
}

template <typename TAlgebra>
double NewtonSolver<TAlgebra>::average_linear_steps(size_t call) const
{
	return m_vTotalLinSolverSteps[call] / ((double)m_vLinSolverCalls[call]);
}

template <typename TAlgebra>
int NewtonSolver<TAlgebra>::total_linsolver_calls() const
{
	int allCalls = 0;
	for(size_t call = 0; call < m_vLinSolverCalls.size(); ++call)
		allCalls += m_vLinSolverCalls[call];
	return allCalls;
}

template <typename TAlgebra>
int NewtonSolver<TAlgebra>::total_linsolver_steps() const
{
	int allSteps = 0;
	for(size_t call = 0; call < m_vLinSolverCalls.size(); ++call)
		allSteps += m_vTotalLinSolverSteps[call];
	return allSteps;
}

template <typename TAlgebra>
double NewtonSolver<TAlgebra>::total_average_linear_steps() const
{
	return total_linsolver_steps()/((double)total_linsolver_calls());
}


template <typename TAlgebra>
void NewtonSolver<TAlgebra>::clear_average_convergence()
{
	m_vLinSolverRates.clear();
	m_vNonLinSolverRates.clear();
	m_vLinSolverCalls.clear();
	m_vTotalLinSolverSteps.clear();
}

template <typename TAlgebra>
void NewtonSolver<TAlgebra>::write_debug(const vector_type& vec, const char* filename)
{
//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
	name.append(ext).append(".vec");

//	write
	typedef DebugWritingObject<TAlgebra> base_writer_type;
	base_writer_type::write_debug(vec, name.c_str());
}

template <typename TAlgebra>
void NewtonSolver<TAlgebra>::write_debug(const matrix_type& mat, const char* filename)
{
//	add iter count to name
	std::string name(filename);
	char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
	name.append(ext).append(".mat");

//	write
	typedef DebugWritingObject<TAlgebra> base_writer_type;
	base_writer_type::write_debug(mat, name.c_str());
}

template <typename TAlgebra>
std::string NewtonSolver<TAlgebra>::config_string() const
{
	std::stringstream ss;
	ss << "NewtonSolver\n";
	ss << " LinearSolver: ";
	if(m_spLinearSolver.valid())	ss << ConfigShift(m_spLinearSolver->config_string()) << "\n";
	else							ss << " NOT SET!\n";
	ss << " ConvergenceCheck: ";
	if(m_spConvCheck.valid())		ss << ConfigShift(m_spConvCheck->config_string()) << "\n";
	else							ss << " NOT SET!\n";
	ss << " LineSearch: ";
	if(m_spLineSearch.valid())		ss << ConfigShift(m_spLineSearch->config_string()) << "\n";
	else							ss << " not set.\n";
	return ss.str();
}


}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__ */


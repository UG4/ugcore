/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel, modifications nested Newton: Markus Knodel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON_IMPL__

#include "newton.h"

#include <iostream>
#include <sstream>
// #include <limits>

// #include "lib_disc/function_spaces/grid_function_util.h"
#include "common/util/string_util.h"

//#include "lib_disc/operator/non_linear_operator/newton_solver/nestedNewtonRFSwitch.h"

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

namespace ug {

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver,
			SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
			SmartPtr<ILineSearch<vector_type> > spLineSearch) :
			m_spLinearSolver(LinearSolver),
			m_spConvCheck(spConvCheck),
			m_spLineSearch(spLineSearch),
			m_N(nullptr),
			m_J(nullptr),
			m_spAss(nullptr),
			m_reassembe_J_freq(0),
			m_dgbCall(0),
			m_lastNumSteps(0),
			m_newtonUpdater(nullptr)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//			m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
{};

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver() :
	m_spLinearSolver(nullptr),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_spLineSearch(nullptr),
	m_N(nullptr),
	m_J(nullptr),
	m_spAss(nullptr),
	m_reassembe_J_freq(0),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_newtonUpdater(nullptr)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//	,
//	m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
{};

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver(SmartPtr<IOperator<vector_type> > N) :
	m_spLinearSolver(nullptr),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_spLineSearch(nullptr),
	m_N(nullptr),
	m_J(nullptr),
	m_spAss(nullptr),
	m_reassembe_J_freq(0),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_newtonUpdater(nullptr)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//	,
//	m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
{
	init(N);
};

template <typename TAlgebra>
NewtonSolver<TAlgebra>::
NewtonSolver(SmartPtr<IAssemble<TAlgebra> > spAss) :
	m_spLinearSolver(nullptr),
	m_spConvCheck(new StdConvCheck<vector_type>(10, 1e-8, 1e-10, true)),
	m_spLineSearch(nullptr),
	m_N(nullptr),
	m_J(nullptr),
	m_spAss(nullptr),
	m_reassembe_J_freq(0),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_newtonUpdater(nullptr)
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//	,
//	m_newtonUpdater(new NewtonUpdaterGeneric<vector_type>{})
//#endif
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

//	loop counts (for the the convergence rate statistics etc.)
	int loopCnt = 0;
	m_lastNumSteps = 0;
	
//	write start defect for debug
	char debug_name_ext[20];
	if (this->debug_writer_valid())
	{
		snprintf(debug_name_ext, 20, "_iter%03d", loopCnt);
		write_debug(*spD, std::string("NEWTON_Defect") + debug_name_ext);
		write_debug(u, "NEWTON_StartSolution");
	}

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
		m_lastNumSteps = loopCnt;

	// 	set c = 0
		NEWTON_PROFILE_BEGIN(NewtonSetCorretionZero);
		spC->set(0.0);
		NEWTON_PROFILE_END();

		for(size_t i = 0; i < m_innerStepUpdate.size(); ++i)
			m_innerStepUpdate[i]->update();

	// 	Compute Jacobian
		try{
			if(m_reassembe_J_freq == 0 || loopCnt % m_reassembe_J_freq == 0) // if we need to reassemble
			{
				NEWTON_PROFILE_BEGIN(NewtonComputeJacobian);
				m_J->init(u);
				NEWTON_PROFILE_END();
			}
		}UG_CATCH_THROW("NewtonSolver::apply: Initialization of Jacobian failed.");

	//	Write the current Jacobian for debug and prepare the section for the lin. solver
		if (this->debug_writer_valid())
		{
			write_debug(m_J->get_matrix(), std::string("NEWTON_Jacobian") + debug_name_ext);
			this->enter_debug_writer_section(std::string("NEWTON_LinSolver") + debug_name_ext);
		}

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

		this->leave_debug_writer_section();
		
	//	store convergence history
		const int numSteps = m_spLinearSolver->step();
		if(loopCnt >= static_cast<int>(m_vTotalLinSolverSteps.size())) m_vTotalLinSolverSteps.resize(loopCnt+1);
		if(loopCnt >= static_cast<int>(m_vLinSolverCalls.size())) m_vLinSolverCalls.resize(loopCnt+1, 0);
		if(loopCnt >= static_cast<int>(m_vLinSolverRates.size())) m_vLinSolverRates.resize(loopCnt+1, 0);
		m_vTotalLinSolverSteps[loopCnt] += numSteps;
		m_vLinSolverCalls[loopCnt] += 1;
		m_vLinSolverRates[loopCnt] += m_spLinearSolver->convergence_check()->avg_rate();

		try{
		// 	Line Search
			if(m_spLineSearch.valid())
			{
				m_spLineSearch->set_offset("   #  ");
				NEWTON_PROFILE_BEGIN(NewtonLineSearch);

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE

				if( m_newtonUpdater != nullptr )
					m_spLineSearch->setNewtonUpdater(m_newtonUpdater);

//#endif
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
				if( m_newtonUpdater != nullptr )
				{
//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE

					if( ! m_newtonUpdater->updateSolution(u,*spC) )
					{
						UG_LOG("ERROR in 'NewtonSolver::apply': "
								"Newton Update did not work.\n");
						// TODO FIXME was macht conv check update? wie kriege ich hier einfach ein riesiges Residuum rein, ohne was zu rechnen?
						return false;
					}

					if( ! m_newtonUpdater->tellAndFixUpdateEvents(u) )
					{
						UG_LOG("unable to fix local Newton updates" << std::endl );
						return false;
					}

				}
				else
				{
//#else
					// 	update solution
					u -= *spC;
				}
//#endif

			// 	compute new Defect
				NEWTON_PROFILE_BEGIN(NewtonComputeDefect);
				m_N->prepare(u);
				m_N->apply(*spD, u);
				NEWTON_PROFILE_END();
			}
		}UG_CATCH_THROW("NewtonSolver::apply: Line Search update failed.");


	//	update counter
		loopCnt++;

	// 	check convergence
		m_spConvCheck->update(*spD);
		if(loopCnt-1 >= static_cast<int>(m_vNonLinSolverRates.size())) m_vNonLinSolverRates.resize(loopCnt, 0);
		m_vNonLinSolverRates[loopCnt-1] += m_spConvCheck->rate();

	//	write defect for debug
		if (this->debug_writer_valid())
		{
			snprintf(debug_name_ext, 20, "_iter%03d", loopCnt);
			write_debug(*spD, std::string("NEWTON_Defect") + debug_name_ext);
			write_debug(*spC, std::string("NEWTON_Correction") + debug_name_ext);
			write_debug(u, std::string("NEWTON_Solution") + debug_name_ext);
		}
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
		UG_LOG(std::setw(13) << std::setprecision(2) << std::fixed << m_vTotalLinSolverSteps[call] / static_cast<double>(m_vLinSolverCalls[call]) << " | ");
		allNonLinRatesProduct *= pow(static_cast<number>(m_vNonLinSolverRates[call])/static_cast<double>(m_vLinSolverCalls[call]),static_cast<double>(m_vLinSolverCalls[call]));
		UG_LOG(std::setw(16) << std::setprecision(6) << std::scientific << m_vNonLinSolverRates[call] / static_cast<double>(m_vLinSolverCalls[call]) << " | ");
		allLinRatesProduct *= (number)std::pow(static_cast<number>(m_vLinSolverRates[call])/static_cast<double>(m_vLinSolverCalls[call]),static_cast<number>(m_vTotalLinSolverSteps[call]));
		UG_LOG(std::setw(13) << std::setprecision(6) << std::scientific << m_vLinSolverRates[call] / static_cast<double>(m_vLinSolverCalls[call]));
		UG_LOG("\n");
	}
	UG_LOG( "        all | ");
	UG_LOG(std::setw(9) << allCalls << " | ");
	UG_LOG(std::setw(15) << allLinSteps << " | ");
	UG_LOG(std::setw(13) << std::setprecision(2) << std::fixed << allLinSteps / static_cast<number>(allCalls) << " | ");
	UG_LOG(std::setw(16) << std::setprecision(6) << std::scientific << std::pow((number)allNonLinRatesProduct,(number)1.0/static_cast<number>(allCalls)) << " | ");
	UG_LOG(std::setw(13) << std::setprecision(6) << std::scientific << std::pow((number)allLinRatesProduct,(number)1.0/static_cast<number>(allLinSteps)));
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
	return m_vTotalLinSolverSteps[call] / static_cast<double>(m_vLinSolverCalls[call]);
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
	return total_linsolver_steps()/static_cast<double>(total_linsolver_calls());
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
void NewtonSolver<TAlgebra>::write_debug(const vector_type& vec, std::string name)
{
//	add call count to name
	char ext[20]; snprintf(ext, 20, "_call%03d", m_dgbCall);

//	write
	using base_writer_type = DebugWritingObject<TAlgebra>;
	base_writer_type::write_debug(vec, name + ext);
}

template <typename TAlgebra>
void NewtonSolver<TAlgebra>::write_debug(const matrix_type& mat, std::string name)
{
//	add call count to name
	char ext[20]; snprintf(ext, 20, "_call%03d", m_dgbCall);

//	write
	using base_writer_type = DebugWritingObject<TAlgebra>;
	base_writer_type::write_debug(mat, name + ext);
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
	if(m_reassembe_J_freq != 0)		ss << " Reassembling Jacobian only once per " << m_reassembe_J_freq << " step(s)\n";
	return ss.str();
}


}

#endif
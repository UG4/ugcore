/*
 * Copyright (c) 2010 - 2017:  G-CSC, Goethe University Frankfurt
 * Author: Arne Naegel
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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NI_SOLVER__NI_IMPL__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NI_SOLVER__NI_IMPL__

#include <iostream>
#include <sstream>

#include "nested_iteration.h"
#include "lib_disc/function_spaces/grid_function_util.h"
#include "common/util/string_util.h"

#define PROFILE_NESTED_ITER
#ifdef PROFILE_NESTED_ITER
	#define NESTED_ITER_PROFILE_FUNC()		PROFILE_FUNC_GROUP("NestedIteration")
	#define NESTED_ITER_PROFILE_BEGIN(name)	PROFILE_BEGIN_GROUP(name, "NestedIteration")
	#define NESTED_ITER_PROFILE_END()		PROFILE_END()
#else
	#define NESTED_ITER_PROFILE_FUNC()
	#define NESTED_ITER_PROFILE_BEGIN(name)
	#define NESTED_ITER_PROFILE_END()
#endif

namespace ug{

template <typename TDomain, typename TAlgebra>
NestedIterationSolver<TDomain,TAlgebra>::
NestedIterationSolver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver) :
			m_spLinearSolver(LinearSolver),
			m_N(NULL),
			m_J(NULL),
			m_spAss(NULL),
			m_dgbCall(0),
			m_lastNumSteps(0),
			m_bUseAdaptiveRefinement(true),
			m_maxSteps(100), m_TOL(1e-3), m_absTOL(1e-12)
{};

template <typename TDomain, typename TAlgebra>
NestedIterationSolver<TDomain,TAlgebra>::
NestedIterationSolver() :
	m_spLinearSolver(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(NULL),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_bUseAdaptiveRefinement(true),
	m_maxSteps(100), m_TOL(1e-3), m_absTOL(1e-12)
{};

template <typename TDomain, typename TAlgebra>
NestedIterationSolver<TDomain,TAlgebra>::
NestedIterationSolver(SmartPtr<IOperator<vector_type> > N) :
	m_spLinearSolver(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(NULL),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_bUseAdaptiveRefinement(true),
	m_maxSteps(100), m_TOL(1e-3), m_absTOL(1e-12)
{
	init(N);
};

template <typename TDomain, typename TAlgebra>
NestedIterationSolver<TDomain,TAlgebra>::
NestedIterationSolver(SmartPtr<IAssemble<TAlgebra> > spAss) :
	m_spLinearSolver(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(spAss),
	m_spDomErr(spAss),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_baseLevel(0),
	m_bUseAdaptiveRefinement(true),
	m_maxSteps(100), m_TOL(1e-3), m_absTOL(1e-12)
{
	m_N = SmartPtr<AssembledOperator<TAlgebra> >(new AssembledOperator<TAlgebra>(m_spAss));
};


template <typename TDomain, typename TAlgebra>
NestedIterationSolver<TDomain,TAlgebra>::
NestedIterationSolver(SmartPtr<IAssemble<TAlgebra> > spAss, SmartPtr<IAssemble<TAlgebra> > spDomErr) :
	m_spLinearSolver(NULL),
	m_N(NULL),
	m_J(NULL),
	m_spAss(spAss),
	m_spDomErr(spDomErr),
	m_dgbCall(0),
	m_lastNumSteps(0),
	m_baseLevel(0),
	m_bUseAdaptiveRefinement(true)
{
	m_N = SmartPtr<AssembledOperator<TAlgebra> >(new AssembledOperator<TAlgebra>(m_spAss));
};


/*! Requires an AssembledOperator */
template <typename TDomain, typename TAlgebra>
bool NestedIterationSolver<TDomain,TAlgebra>::init(SmartPtr<IOperator<vector_type> > N)
{
	NESTED_ITER_PROFILE_BEGIN(NestedIterationSolver_init);
	m_N = N.template cast_dynamic<AssembledOperator<TAlgebra> >();
	if(m_N.invalid())
		UG_THROW("NestedIterationSolver: currently only works for AssembledDiscreteOperator.");

	m_spAss = m_N->discretization();
	return true;
}

//!	todo: remove from interface
template <typename TDomain, typename TAlgebra>
bool NestedIterationSolver<TDomain,TAlgebra>::prepare(vector_type& u)
{

	return true;
}


//! Refines domain and provides error estimate.
/*! Values depend on m_spDomErr */
template <typename TDomain, typename TAlgebra>
void NestedIterationSolver<TDomain,TAlgebra>::estimate_and_mark_domain(const grid_function_type& u, SmartPtr<IElementMarkingStrategy<TDomain> > spMarking, bool bClearMarks)
{
	  UG_LOG("Computing error... "<< std::endl);

	  typedef DomainDiscretization<TDomain, TAlgebra> domain_disc_type;
	 //  typedef IDomainErrorIndicator<TAlgebra> domain_indicator_type;

	  SmartPtr<domain_disc_type> spDomainEstimator = m_spDomErr.template cast_dynamic<domain_disc_type>();

	  // some checks
	  UG_ASSERT(spDomainEstimator.valid(), "Huhh: No estimatable object???")
	  UG_ASSERT(m_spRefiner.valid(), "Huhh: Invalid refiner???");

	  // compute error
	  if (m_spElemError.valid()) {
		  // debug version
		  spDomainEstimator->calc_error(u, u.dd(), m_spElemError.get());
	  } else {
		  // standard version
		  spDomainEstimator->calc_error(u, u.dd());
	  }

	  // set (new) marks
	  if (bClearMarks) m_spRefiner->clear_marks();
	  spDomainEstimator->mark_with_strategy(*m_spRefiner, spMarking);
	  UG_LOG("Estimated error: " << spMarking->global_estimated_error());




}


//! Coarsen all elements
template <typename TDomain, typename TAlgebra>
number NestedIterationSolver<TDomain,TAlgebra>::coarsen_domain(const grid_function_type& u)
{
	/*typedef typename domain_traits<TDomain::dim>::element_type TElem;
	SmartPtr<DoFDistribution> spDD=u.dof_distribution();
	m_spRefiner->mark(spDD->begin<TElem>(), spDD->end<TElem>(), RM_COARSEN);
	m_spRefiner->coarsen();*/
	
	return 0; // dummy value
}


//! Apply solver for top level grid function
/*! returns: approximation by nested iteration. Input values are ignored.
 *
 * */
template <typename TDomain, typename TAlgebra>
bool NestedIterationSolver<TDomain,TAlgebra>::apply(vector_type& u)
{
#ifdef UG_PARALLEL
	// UG_ASSERT(0, "Not implemented!")
#endif
	NESTED_ITER_PROFILE_BEGIN(NestedIterationSolver_apply);

	//	increase call count
	m_dgbCall++;
	grid_function_type &usol = dynamic_cast<grid_function_type&>(u);

	UG_ASSERT (usol.grid_level().is_surface(), "Must supply surface grid function");
	UG_ASSERT (usol.grid_level().top(), "Must supply top level grid function");
	const GridLevel& surfGridLevel = usol.grid_level();

	// Check for linear solver.
	if(m_spLinearSolver.invalid())
		UG_THROW("NestedIterationSolver::apply: Linear Solver not set.");

	// Check for Jacobian.
	if(m_J.invalid() || m_J->discretization() != m_spAss) {
		m_J = make_sp(new AssembledLinearOperator<TAlgebra>(m_spAss));
	}
	m_J->set_level(m_N->level());

	UG_LOG("N level: "<<m_N->level()<< std::endl);
	UG_LOG("u level: "<< usol.grid_level() << std::endl);

	// Create tmp vectors
	SmartPtr<vector_type> spD = u.clone_without_values();

	char ext[50];
	int loopCnt = 0;
	m_lastNumSteps = 0;
	{
		//	write start defect for debug
		sprintf(ext, "_call%03d", m_dgbCall);
		std::string name("NESTED_ITER_StartSolution");
		name.append(ext);
		write_debug(u, name.c_str());
	}

	// Assembly loop (from some coarse level to TOP)
	int surfLevel = usol.grid_level().level();  // todo: should start on coarse(r) levels
	m_topLevel = surfLevel+m_maxSteps;
	for (int lev=surfLevel; lev<m_topLevel; ++lev)
	{
		/* OLD LUA CODE:
		 *  globalDisc:assemble_linear(A, b)
		 * 	globalDisc:adjust_solution(u)
		 * 	solver:init(A)
		 * 	solver:apply_return_defect(u,b)
		 * 	// globalDisc:adjust_solution(u)
		 */

		// Assemble.
		NESTED_ITER_PROFILE_BEGIN(NestedIterationAssemble);
		m_spAss->assemble_linear(*m_J, *spD, surfGridLevel);   // todo: replace for non-linear problems
		m_spAss->adjust_solution(u, surfGridLevel);
		NESTED_ITER_PROFILE_END();

		sprintf(ext, "_call%03d_iter%03d", m_dgbCall, loopCnt);
		{
			// write initial data for debug
			std::string name0("NESTED_ITER_InitialSolution"); name0.append(ext);
			write_debug(u, name0.c_str());
			std::string name("NESTED_ITER_InitialDefect"); name.append(ext);
			write_debug(*spD, name.c_str());

			//	Write Jacobian for debug
			std::string matname("NESTED_ITER_Jacobian");
			matname.append(ext);
			write_debug(m_J->get_matrix(), matname.c_str());
		}

		// Initialize inverse of Jacobian.
		try{
			NESTED_ITER_PROFILE_BEGIN(NestedIterationPrepareLinSolver);
			if(!m_spLinearSolver->init(m_J, u))
			{
				UG_LOG("ERROR in 'NestedIterationSolver::apply': Cannot init Inverse Linear "
						"Operator for Jacobi-Operator.\n");
				return false;
			}
			NESTED_ITER_PROFILE_END();
		}UG_CATCH_THROW("NestedIterationSolver::apply: Initialization of Linear Solver failed.");

		// Solve linearized system.
		try{
			NESTED_ITER_PROFILE_BEGIN(NewtonApplyLinSolver);
			if(!m_spLinearSolver->apply(u, *spD))
			{
				UG_LOG("ERROR in 'NestedIterationSolver::apply': Cannot apply Inverse Linear "
						"Operator for Jacobi-Operator.\n");
				return false;
			}
			NESTED_ITER_PROFILE_END();
		}UG_CATCH_THROW("NestedIterationSolver::apply: Application of Linear Solver failed.");

		// Adjust solution.
		m_spAss->adjust_solution(u, surfGridLevel);
		{
			//	write defect for debug
			std::string name("NESTED_ITER_FinalDefect"); name.append(ext);
			write_debug(*spD, name.c_str());
			std::string name3("NESTED_ITER_FinalSolution"); name3.append(ext);
			write_debug(u, name3.c_str());
		}

		// Update counter.
		loopCnt++;

		// Quit?
		if (use_adaptive_refinement() == false) {
			UG_LOG("SUCCESS: Non-adaptive version!" <<  std::endl);
			break;
		}

		// Estimate and mark domain.
		this->estimate_and_mark_domain(usol, m_spRefinementMarking);

		//	OPTIONAL: write defect for debug
		if (m_spElemError.valid())
		{
			std::string name("NESTED_ITER_Error");
			name.append(ext);
			VTKOutput<TDomain::dim> outError;
			outError.print(name.c_str(), *m_spElemError);
		}

		// Check error.
		const number err = m_spRefinementMarking->global_estimated_error();
		double desiredTOL = (m_spAssociatedSpace.valid()) ? m_TOL*m_spAssociatedSpace->norm(usol) + m_absTOL : m_TOL;

		UG_LOG("NI *** Desired tolerance: " << m_TOL << " * "<<  m_spAssociatedSpace->norm(usol) <<  " + "<< m_absTOL <<std::endl);

		// Quit or continue?
		if(err < desiredTOL)
		{
			// Quit.
			UG_LOG("SUCCESS: Error smaller than tolerance: " << err << " < "<< desiredTOL << std::endl);
			break;
		} else {
			// Refine and continue.
			UG_LOG("FAILED: Error (still) greater than tolerance: " << err << " > "<< desiredTOL << std::endl);
			m_spRefiner->refine();
		}
	}


	return true;
}


template <typename TDomain, typename TAlgebra>
std::string NestedIterationSolver<TDomain,TAlgebra>::config_string() const
{
	std::stringstream ss;
	ss << "NestedIterationSolver\n";
	ss << " LinearSolver: ";
	if(m_spLinearSolver.valid())	ss << ConfigShift(m_spLinearSolver->config_string()) << "\n";
	else							ss << " NOT SET!\n";
	return ss.str();
}


}

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NESTED_ITER_SOLVER__NESTED_ITER_IMPL__ */


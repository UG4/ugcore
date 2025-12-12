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

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__

//#include <cmath>

#include "lib_algebra/operator/interface/operator_inverse.h"
#include "lib_algebra/operator/interface/linear_operator_inverse.h"
#include "lib_algebra/operator/debug_writer.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "../line_search.h"
#include "newton_update_interface.h"
// #include "lib_algebra/operator/debug_writer.h"

//#include "nestedNewtonRFSwitch.h"

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
//
//#include "newtonUpdaterGeneric.h"
//
//#endif

namespace ug {

/// Newton solver for assembling based discretizations
template <typename TAlgebra>
class NewtonSolver
	: 	public IOperatorInverse<typename TAlgebra::vector_type>,
		public DebugWritingObject<TAlgebra>
{
	public:
	///	Algebra type
		using algebra_type = TAlgebra;

	///	Vector type
		using vector_type = typename TAlgebra::vector_type;

	///	Matrix type
		using matrix_type = typename TAlgebra::matrix_type;

	public:
	///	constructor setting operator
	explicit NewtonSolver(SmartPtr<IOperator<vector_type> > N);

	///	constructor using assembling
	explicit NewtonSolver(SmartPtr<IAssemble<TAlgebra> > spAss);

	///	default constructor
		NewtonSolver();

	///	constructor
		NewtonSolver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver,
		             SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
		             SmartPtr<ILineSearch<vector_type> > spLineSearch);

		~NewtonSolver() override = default;

	///	sets the linear solver
		void set_linear_solver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver) {m_spLinearSolver = LinearSolver;}

	/// sets the convergence check
		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

	///	sets the line search
		void set_line_search(SmartPtr<ILineSearch<vector_type> > spLineSearch) {m_spLineSearch = spLineSearch;}
		void disable_line_search() {m_spLineSearch = nullptr;}
		SmartPtr<ILineSearch<vector_type> > line_search() {return m_spLineSearch;}

	/// This operator inverts the Operator N: Y -> X
		bool init(SmartPtr<IOperator<vector_type> > N) override;

	/// prepare Operator
		bool prepare(vector_type& u) override;

	/// apply Operator, i.e. N^{-1}(0) = u
		bool apply(vector_type& u) override;

	///	returns information about configuration parameters
		/**
		 * this should return necessary information about parameters and possibly
		 * calling config_string of subcomponents.
		 *
		 * \returns std::string	necessary information about configuration parameters
		 */

		std::string config_string() const override;

	/// prints average linear solver convergence
		void print_average_convergence() const;

	///	information on convergence history
	/// \{
		size_t num_newton_steps() const;
		int num_linsolver_calls(size_t call) const;
		int num_linsolver_steps(size_t call) const;
		double average_linear_steps(size_t call) const;
		int total_linsolver_calls() const;
		int total_linsolver_steps() const;
		double total_average_linear_steps() const;
		int last_num_newton_steps() const	{return m_lastNumSteps;}
	/// \}

	/// resets average linear solver convergence
		void clear_average_convergence();

	///	add inner step update (applied before every linear solver step)
		void add_inner_step_update(SmartPtr<INewtonUpdate > NU)
			{m_innerStepUpdate.push_back(NU);}

	///	clears inner step update
		void clear_inner_step_update(SmartPtr<INewtonUpdate > NU)
			{m_innerStepUpdate.clear();}

	///	add outer step update (applied before every Newton step)
		void add_step_update(SmartPtr<INewtonUpdate > NU)
			{m_stepUpdate.push_back(NU);}

	///	clears outer step update
		void clear_step_update(SmartPtr<INewtonUpdate > NU)
			{m_stepUpdate.clear();}
		
	///	sets the frequency of reassembling of the Jacobian (0 == 1 == in every step, i.e. classically)
		void set_reassemble_J_freq(int freq)
			{m_reassembe_J_freq = freq;};

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE
		void setNewtonUpdater( SmartPtr<NewtonUpdaterGeneric<vector_type> > nU )
		{
			m_newtonUpdater = nU;
		}
//#endif

		bool createNewtonUpdater()
		{
			if( m_newtonUpdater != nullptr )
			{
				m_newtonUpdater = SmartPtr<NewtonUpdaterGeneric<vector_type> >
										  (new NewtonUpdaterGeneric<vector_type>{});

				return true;
			}
			return false;
		}


	private:
	///	help functions for debug output
	///	\{
		void write_debug(const vector_type& vec, std::string filename) override;
		void write_debug(const matrix_type& mat, std::string filename) override;
	/// \}

	private:
	///	linear solver
		SmartPtr<ILinearOperatorInverse<vector_type> > m_spLinearSolver;

	/// Convergence Check
		SmartPtr<IConvergenceCheck<vector_type> > m_spConvCheck;

	/// LineSearch
		SmartPtr<ILineSearch<vector_type> > m_spLineSearch;

	/// Update
		std::vector<SmartPtr<INewtonUpdate> > m_innerStepUpdate;
		std::vector<SmartPtr<INewtonUpdate> > m_stepUpdate;

	///	assembling routine
		SmartPtr<AssembledOperator<algebra_type> > m_N;
	///	jacobi operator
		SmartPtr<AssembledLinearOperator<algebra_type> > m_J;
	///	assembling
		SmartPtr<IAssemble<TAlgebra> > m_spAss;
	/// how often to reassemble the Jacobian (0 == 1 == in every step, i.e. classically)
		int m_reassembe_J_freq;

	///	call counter
		int m_dgbCall;
		int m_lastNumSteps;

	/// convergence history of linear solver
	/// \{
		std::vector<int> m_vTotalLinSolverSteps;
		std::vector<int> m_vLinSolverCalls;
		std::vector<number> m_vNonLinSolverRates;
		std::vector<number> m_vLinSolverRates;
	/// \}

//#if ENABLE_NESTED_NEWTON_RESOLFUNC_UPDATE

		SmartPtr<NewtonUpdaterGeneric<vector_type> > m_newtonUpdater;

//#endif

};

}

#include "newton_impl.h"

#endif
/*
 * newton.h
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__

#include <cmath>

#include "lib_algebra/operator/interface/operator_inverse.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "lib_disc/operator/linear_operator/assembled_linear_operator.h"
#include "../line_search.h"
#include "newton_update_interface.h"
#include "lib_algebra/operator/debug_writer.h"

namespace ug {

/// Newton solver for assembling based discretizations
template <typename TAlgebra>
class NewtonSolver
	: 	public IOperatorInverse<typename TAlgebra::vector_type>,
		public DebugWritingObject<TAlgebra>
{
	public:
	///	Algebra type
		typedef TAlgebra algebra_type;

	///	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	///	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;
		using base_writer_type::write_debug;

	public:
	///	constructor setting operator
		NewtonSolver(SmartPtr<IOperator<vector_type> > N);

	///	constructor using assembling
		NewtonSolver(SmartPtr<IAssemble<TAlgebra> > spAss);

	///	default constructor
		NewtonSolver();

	///	constructor
		NewtonSolver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver,
		             SmartPtr<IConvergenceCheck<vector_type> > spConvCheck,
		             SmartPtr<ILineSearch<vector_type> > spLineSearch);

	///	sets the linear solver
		void set_linear_solver(SmartPtr<ILinearOperatorInverse<vector_type> > LinearSolver) {m_spLinearSolver = LinearSolver;}

	/// sets the convergence check
		void set_convergence_check(SmartPtr<IConvergenceCheck<vector_type> > spConvCheck);

	///	sets the line search
		void set_line_search(SmartPtr<ILineSearch<vector_type> > spLineSearch) {m_spLineSearch = spLineSearch;}

	/// This operator inverts the Operator N: Y -> X
		virtual bool init(SmartPtr<IOperator<vector_type> > N);

	/// prepare Operator
		virtual bool prepare(vector_type& u);

	/// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

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

	private:
	///	help functions for debug output
	///	\{
		void write_debug(const vector_type& vec, const char* filename);
		void write_debug(const matrix_type& mat, const char* filename);
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

	/// line search parameters
	/// \{
		int m_maxLineSearch;
		number m_lambda_start;
		number m_lambda_reduce;
	/// \}

	///	call counter
		int m_dgbCall;

	/// convergence history of linear solver
	/// \{
		std::vector<int> m_vTotalLinSolverSteps;
		std::vector<int> m_vLinSolverCalls;
		std::vector<number> m_vNonLinSolverRates;
		std::vector<number> m_vLinSolverRates;
	/// \}
};

}

#include "newton_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__ */

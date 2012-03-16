/*
 * newton.h
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__
#define __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__

#include <cmath>

#include "lib_algebra/operator/operator_inverse_interface.h"

// modul intern headers
#include "lib_disc/assemble_interface.h"
#include "lib_disc/operator/non_linear_operator/assembled_non_linear_operator.h"
#include "../line_search.h"
#include "lib_algebra/operator/debug_writer.h"

namespace ug {

template <typename TAlgebra>
class NewtonSolver
	: public IOperatorInverse<	typename TAlgebra::vector_type,
	  	  	  	  	  	  	  	typename TAlgebra::vector_type>,
	public DebugWritingObject<TAlgebra>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	Matrix type
		typedef typename TAlgebra::matrix_type matrix_type;

	protected:
		typedef DebugWritingObject<TAlgebra> base_writer_type;
		using base_writer_type::write_debug;

	public:
		NewtonSolver(ILinearOperatorInverse<vector_type, vector_type>& LinearSolver,
					SmartPtr<IConvergenceCheck> spConvCheck,
					SmartPtr<ILineSearch<vector_type> > spLineSearch, bool reallocate) :
					m_pLinearSolver(&LinearSolver),
					m_spConvCheck(spConvCheck),
					m_spLineSearch(spLineSearch),
					m_reallocate(reallocate), m_allocated(false),
					m_dgbCall(0)
			{};

		NewtonSolver() :
			m_pLinearSolver(NULL),
			m_spConvCheck(new StandardConvCheck(10, 1e-8, 1e-10, true)),
			m_spLineSearch(NULL),
			m_reallocate(false), m_allocated(false),
			m_dgbCall(0)
			{};

		void set_linear_solver(ILinearOperatorInverse<vector_type, vector_type>& LinearSolver) {m_pLinearSolver = &LinearSolver;}
		void set_convergence_check(SmartPtr<IConvergenceCheck> spConvCheck)
		{
			m_spConvCheck = spConvCheck;
			m_spConvCheck->set_offset(3);
			m_spConvCheck->set_symbol('#');
			m_spConvCheck->set_name("Newton Solver");
		}

		void set_line_search(SmartPtr<ILineSearch<vector_type> > spLineSearch) {m_spLineSearch = spLineSearch;}

		// init: This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<vector_type, vector_type>& N);

		// prepare Operator
		virtual bool prepare(vector_type& u);

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

		~NewtonSolver();

		// prints average linear solver convergence
		void print_average_convergence() const;

		// resets average linear solver convergence
		void clear_average_convergence();

	private:
		bool allocate_memory(const vector_type& u);
		bool deallocate_memory();

		void write_debug(const vector_type& vec, const char* filename)
		{
		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
			name.append(ext).append(".vec");

		//	write
			base_writer_type::write_debug(vec, name.c_str());
		}

		void write_debug(const matrix_type& mat, const char* filename)
		{
		//	add iter count to name
			std::string name(filename);
			char ext[20]; sprintf(ext, "_call%03d", m_dgbCall);
			name.append(ext).append(".mat");

		//	write
			base_writer_type::write_debug(mat, name.c_str());
		}

	private:
		ILinearOperatorInverse<vector_type, vector_type>* m_pLinearSolver;

		// Convergence Check
		SmartPtr<IConvergenceCheck> m_spConvCheck;

		// LineSearch
		SmartPtr<ILineSearch<vector_type> > m_spLineSearch;

		vector_type m_d;
		vector_type m_c;

		AssembledOperator<algebra_type>* m_N;
		AssembledLinearOperator<algebra_type>* m_J;
		IAssemble<algebra_type>* m_pAss;

		// line search parameters
		int m_maxLineSearch;
		number m_lambda_start;
		number m_lambda_reduce;

		bool m_reallocate;
		bool m_allocated;

		int m_dgbCall;

		// convergence history of linear solver
		std::vector<int> m_vTotalLinSolverSteps;
		std::vector<int> m_vLinSolverCalls;
};

}

#include "newton_impl.h"

#endif /* __H__UG__LIB_DISC__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__ */

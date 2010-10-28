/*
 * newton.h
 *
 *  Created on: 26.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__
#define __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__

#include <cmath>

// modul intern headers
#include "lib_discretization/assemble.h"
#include "lib_discretization/operator/operator.h"
#include "../line_search.h"

namespace ug {

template <typename TDoFDistribution, typename TAlgebra>
class NewtonSolver : public IOperatorInverse<	typename TAlgebra::vector_type,
												typename TAlgebra::vector_type>
{
	public:
	//	Algebra type
		typedef TAlgebra algebra_type;

	//	Vector type
		typedef typename TAlgebra::vector_type vector_type;

	//	DoFDistribution Type
		typedef TDoFDistribution dof_distribution_type;

	public:
		NewtonSolver(ILinearOperatorInverse<vector_type, vector_type>& LinearSolver,
					IConvergenceCheck& ConvCheck,
					ILineSearch<vector_type>* LineSearch, bool reallocate) :
					m_pLinearSolver(&LinearSolver),
					m_pConvCheck(&ConvCheck),
					m_pLineSearch(LineSearch),
					m_reallocate(reallocate), m_allocated(false)
			{};

		NewtonSolver() :
			m_pLinearSolver(NULL), m_pConvCheck(NULL), m_pLineSearch(NULL), m_reallocate(false), m_allocated(false)
			{};

		void set_linear_solver(ILinearOperatorInverse<vector_type, vector_type>& LinearSolver) {m_pLinearSolver = &LinearSolver;}
		void set_convergence_check(IConvergenceCheck& ConvCheck)
		{
			m_pConvCheck = &ConvCheck;
			m_pConvCheck->set_offset(3);
			m_pConvCheck->set_symbol('#');
			m_pConvCheck->set_name("Newton Solver");
		}
		void set_line_search(ILineSearch<vector_type>& LineSearch) {m_pLineSearch = &LineSearch;}

		// init: This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<vector_type, vector_type>& N);

		// prepare Operator
		virtual bool prepare(vector_type& u);

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(vector_type& u);

		~NewtonSolver();

	private:
		bool allocate_memory(const vector_type& u);
		bool deallocate_memory();

	private:
		ILinearOperatorInverse<vector_type, vector_type>* m_pLinearSolver;

		// Convergence Check
		IConvergenceCheck* m_pConvCheck;

		// LineSearch
		ILineSearch<vector_type>* m_pLineSearch;

		vector_type m_d;
		vector_type m_c;

		AssembledOperator<dof_distribution_type, algebra_type>* m_N;
		AssembledLinearOperator<dof_distribution_type, algebra_type>* m_J;
		IAssemble<dof_distribution_type, algebra_type>* m_pAss;

		// line search parameters
		int m_maxLineSearch;
		number m_lambda_start;
		number m_lambda_reduce;

		bool m_reallocate;
		bool m_allocated;
};

}

#include "newton_impl.h"

#endif /* __H__LIBDISCRETIZATION__OPERATOR__NON_LINEAR_OPERATOR__NEWTON_SOLVER__NEWTON__ */

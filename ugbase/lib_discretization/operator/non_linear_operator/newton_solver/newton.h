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

template <typename TFunction>
class NewtonSolver : public IOperatorInverse<TFunction, TFunction>{
	public:
		typedef TFunction function_type;

	public:
		NewtonSolver(ILinearizedOperatorInverse<function_type, function_type>& LinearSolver,
					ConvergenceCheck<TFunction>& ConvCheck,
					LineSearch<TFunction>* LineSearch,	bool reallocate) :
					m_LinearSolver(LinearSolver),
					m_ConvCheck(ConvCheck),
					m_LineSearch(LineSearch),
					m_reallocate(reallocate), m_allocated(false)
			{};

		// init: This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<function_type, function_type>& N);

		// prepare Operator
		virtual bool prepare(function_type& u);

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(function_type& u);

		~NewtonSolver();

	private:
		bool allocate_memory(const function_type& u);
		bool deallocate_memory();

	private:
		ILinearizedOperatorInverse<function_type, function_type>& m_LinearSolver;

		// Convergence Check
		ConvergenceCheck<TFunction>& m_ConvCheck;

		// LineSearch
		LineSearch<TFunction>* m_LineSearch;

		function_type* m_d;
		function_type* m_c;

		AssembledOperator<function_type>* m_N;
		AssembledLinearizedOperator<function_type>* m_J;
		IAssemble<function_type>* m_ass;

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

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

namespace ug {

template <typename TDiscreteFunction>
class NewtonSolver : public IOperatorInverse<TDiscreteFunction, TDiscreteFunction>{
	public:
		typedef TDiscreteFunction discrete_function_type;

	public:
		NewtonSolver(ILinearizedOperatorInverse<discrete_function_type, discrete_function_type>& LinearSolver,
						int MaxIterations, number absTol, number relTol,
						int MaxLineSearch, number lambda_start, number lambda_reduce, bool reallocate) :
					m_LinearSolver(LinearSolver),
					m_MaxIterations(MaxIterations), m_absTol(absTol), m_relTol(relTol),
					m_maxLineSearch(MaxLineSearch), m_lambda_start(lambda_start), m_lambda_reduce(lambda_reduce),
					m_reallocate(reallocate), m_allocated(false)
			{};

		// init: This operator inverts the Operator N: Y -> X
		virtual bool init(IOperator<discrete_function_type, discrete_function_type>& N);

		// prepare Operator
		virtual bool prepare(discrete_function_type& u);

		// apply Operator, i.e. N^{-1}(0) = u
		virtual bool apply(discrete_function_type& u);

		~NewtonSolver();

	private:
		bool allocate_memory(const discrete_function_type& u);
		bool deallocate_memory();

		bool is_valid_number(number value)
		{
			// (value >= std::numeric_limits<number>::min() ) == true if value > -infty
			// (value <= std::numeric_limits<number>::max() ) == true if value < infty
			// (value == value                         ) == true if value != NaN
			if(value == 0.0) return true;
			else return value >= std::numeric_limits<number>::min() && value <= std::numeric_limits<number>::max() && value == value && value >= 0;
		}

		bool VecScaleAdd(discrete_function_type& a_func, discrete_function_type& b_func, number s)
		{
			typename discrete_function_type::vector_type& a = a_func.get_vector();
			typename discrete_function_type::vector_type& b = b_func.get_vector();
			typename discrete_function_type::algebra_type::matrix_type::local_matrix_type locMat(1, 1);
			typename discrete_function_type::algebra_type::matrix_type::local_index_type locInd(1);
			typename discrete_function_type::vector_type::local_vector_type locVec(1);

			for(size_t i = 0; i < a.size(); ++i){
				locInd[0][0] = i;
				b.get(locVec, locInd);
				locVec[0] *= s;

				a.add(locVec, locInd);
			}
			return true;
		}

	private:
		ILinearizedOperatorInverse<discrete_function_type, discrete_function_type>& m_LinearSolver;
		int m_MaxIterations;
		number m_absTol;
		number m_relTol;
		discrete_function_type* m_d;
		discrete_function_type* m_c;

		AssembledOperator<discrete_function_type>* m_N;
		AssembledLinearizedOperator<discrete_function_type>* m_J;
		IAssemble<discrete_function_type>* m_ass;

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

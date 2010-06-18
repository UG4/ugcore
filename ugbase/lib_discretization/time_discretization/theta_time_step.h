/*
 * theta_time_step.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP__
#define __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// modul intern libraries
#include "lib_discretization/assemble.h"
#include "lib_discretization/time_discretization/time_discretization_interface.h"

namespace ug{

template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class ThetaTimeDiscretization : public ITimeDiscretization<TDiscreteFunction, TAlgebra> {
	public:
		typedef TDiscreteFunction discrete_function_type;

		// type of algebra
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;


	public:
		// theta = 0 -> Backward Euler
		ThetaTimeDiscretization(IDomainDiscretization<discrete_function_type, algebra_type>& sd, number theta);

		// implements the time step interface
		bool prepare_step(std::deque<discrete_function_type*>& u_old, std::deque<number>& time_old, number dt);

		// implements the assemble interface
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u);
		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u);
		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u);
		IAssembleReturn assemble_solution(discrete_function_type& u);
		IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, const discrete_function_type& u);

	private:
		size_t m_previousSteps;
		number s_a[2];
		number s_m[2];
		std::deque<discrete_function_type*> *m_u_old;
		std::deque<number> *m_time_old;
		number m_dt; // current time step
		number m_time_future;
};

}


#include "theta_time_step_impl.h"

#endif /* __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP__ */

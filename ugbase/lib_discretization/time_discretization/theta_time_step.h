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

template <	typename TDoFDistribution,
			typename TAlgebra>
class ThetaTimeDiscretization : public ITimeDiscretization<TDoFDistribution, TAlgebra>
{
	public:
	//	DoF Distribution Type
		typedef TDoFDistribution dof_distribution_type;

	// 	Type of algebra
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
		// theta = 0 -> Backward Euler
		ThetaTimeDiscretization(IDomainDiscretization<dof_distribution_type, algebra_type>& sd, number theta);

		ThetaTimeDiscretization()
		{
			set_theta(1.0);
		}

		void set_theta(number theta)
		{
			s_a[0] = 1.-theta;
			s_a[1] = theta;
			s_m[0] = 1.;
			s_m[1] = -1.;
			m_previousSteps = 1;
		}
		void set_domain_discretization(IDomainDiscretization<dof_distribution_type, algebra_type>& dd) {this->m_dd = &dd;}

		// return number of previous time steps needed
		size_t num_prev_steps() {return m_previousSteps;}

		// implements the time step interface
		bool prepare_step(std::deque<vector_type*>& u_old, std::deque<number>& time_old, number dt);

	public:
		// implements the assemble interface
		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr);
		IAssembleReturn assemble_jacobian(matrix_type& J, const vector_type& u, const dof_distribution_type& dofDistr);
		IAssembleReturn assemble_defect(vector_type& d, const vector_type& u, const dof_distribution_type& dofDistr);
		IAssembleReturn assemble_solution(vector_type& u, const dof_distribution_type& dofDistr);
		IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, const vector_type& u, const dof_distribution_type& dofDistr);

	private:
		size_t m_previousSteps;
		number s_a[2];
		number s_m[2];
		std::deque<vector_type*> *m_u_old;
		std::deque<number> *m_time_old;
		number m_dt; // current time step
		number m_time_future;
};

}


#include "theta_time_step_impl.h"

#endif /* __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP__ */

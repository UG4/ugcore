/*
 * theta_time_step_impl.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP_IMPL__
#define __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP_IMPL__

#include "theta_time_step.h"

namespace ug{

template <	typename TDiscreteFunction, typename TAlgebra >
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
ThetaTimeDiscretization(IDomainDiscretization<discrete_function_type, algebra_type>& sd, number theta)
	: ITimeDiscretization<TDiscreteFunction, TAlgebra>(sd)
{
	s_a[0] = 1.-theta;
	s_a[1] = theta;
	s_m[0] = 1.;
	s_m[1] = -1.;
	m_previousSteps = 1;
}

template <	typename TDiscreteFunction, typename TAlgebra >
bool
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
prepare_step(std::deque<discrete_function_type*>& u_old, std::deque<number>& time_old, number dt)
{
	if(u_old.size() != m_previousSteps)
	{
		UG_LOG("ERROR in prepare_step: Number of previous solutions is not adequate for this time solver" << std::endl);
		return false;
	}
	if(time_old.size() != m_previousSteps)
	{
		UG_LOG("ERROR in prepare_step: Number of previous time steps is not adequate for this time solver" << std::endl);
		return false;
	}
	if(dt < 0.0)
	{
		UG_LOG("ERROR in prepare_step: Time step size can not be negative." << std::endl);
		return false;
	}

	m_u_old = &u_old;
	m_time_old = &time_old;
	m_dt = dt;
	m_time_future = m_dt + (*m_time_old)[0];

	return true;
}

template <	typename TDiscreteFunction, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u)
{
	// future solution part
	if(this->m_dd.assemble_defect(d, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;
	if(this->m_dd.assemble_jacobian(J, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;

	// previous time step part
	for(size_t i=0; i < m_previousSteps; ++i)
	{
		if(this->m_dd.assemble_defect(d, *(*m_u_old)[i], (*m_time_old)[i], s_m[i+1], s_a[i+1]*m_dt) == IAssemble_OK) return IAssemble_ERROR;
	}

	return IAssemble_OK;
}

template <	typename TDiscreteFunction, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
assemble_jacobian(matrix_type& J, const discrete_function_type& u)
{
	if(this->m_dd.assemble_jacobian(J, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;

	return IAssemble_OK;
}

template <	typename TDiscreteFunction, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
assemble_defect(vector_type& d, const discrete_function_type& u)
{
	// future solution part
	if(this->m_dd.assemble_defect(d, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;

	// previous time step part
	for(size_t i=0; i < m_previousSteps; ++i)
	{
		if(this->m_dd.assemble_defect(d, *(*m_u_old)[i], (*m_time_old)[i], s_m[i+1], s_a[i+1]*m_dt) != IAssemble_OK) return IAssemble_ERROR;
	}

	return IAssemble_OK;
}

template <	typename TDiscreteFunction, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
assemble_solution(discrete_function_type& u)
{
	IAssembleReturn res;

	res = this->m_dd.assemble_solution(u, m_time_future);

	switch(res)
	{
	case IAssemble_ERROR:
			UG_LOG("ERROR in assemble_solution" << std::endl);
			return IAssemble_ERROR;
	case IAssemble_NOT_IMPLEMENTED:
			UG_LOG("ERROR in assemble_solution: function not implemented" << std::endl);
			return IAssemble_ERROR;
	case IAssemble_TIME_INDEPENDENT:
			UG_LOG("ERROR in assemble_solution: Problem time independent" << std::endl);
			return IAssemble_ERROR;

	case IAssemble_OK:
	case IAssemble_NON_LINEAR:	return IAssemble_OK;
	}
	return IAssemble_OK;
}

template <	typename TDiscreteFunction, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDiscreteFunction, TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, const discrete_function_type& u)
{
	return IAssemble_NOT_IMPLEMENTED;
}

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__THETA_TIME_STEP_IMPL__ */

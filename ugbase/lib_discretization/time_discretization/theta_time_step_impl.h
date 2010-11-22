/*
 * theta_time_step_impl.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__THETA_TIME_STEP_IMPL__
#define __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__THETA_TIME_STEP_IMPL__

#include "theta_time_step.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
prepare_step(std::deque<vector_type*>& u_old, std::deque<number>& time_old,
             number dt)
{
//	perform checks
	if(u_old.size() != m_prevSteps)
	{
		UG_LOG("ERROR in ThetaTimeDiscretization::prepare_step:"
			" Number of previous solutions must be "<<m_prevSteps<<".\n");
		return false;
	}
	if(time_old.size() != m_prevSteps)
	{
		UG_LOG("ERROR in ThetaTimeDiscretization::prepare_step: "
			"Number of previous time steps must be "<< m_prevSteps <<".\n");
		return false;
	}
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:prepare_step': "
				"Domain Discretization not set.\n");
		return false;
	}

//	remember old values
	m_pSolOld = &u_old;
	m_pTimeOld = &time_old;

//	remember time step size
	m_dt = dt;

//	compute future time
	m_futureTime = m_dt + (*m_pTimeOld)[0];

//	done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u,
                  const dof_distribution_type& dofDistr)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:assemble_jacobian':"
				" Domain Discretization not set.\n");
		return IAssemble_ERROR;
	}

//	assemble jacobian using current iterate
	if(this->m_pDomDisc->assemble_jacobian
			(J, u, dofDistr, m_futureTime, s_m[0], s_a[0]*m_dt) != IAssemble_OK)
		return IAssemble_ERROR;

//	we're done
	return IAssemble_OK;
}

template <typename TDoFDistribution, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:assemble_defect':"
				" Domain Discretization not set.\n");
		return IAssemble_ERROR;
	}

// 	future solution part
	if(this->m_pDomDisc->assemble_defect
			(d, u, dofDistr, m_futureTime, s_m[0], s_a[0]*m_dt) != IAssemble_OK)
		return IAssemble_ERROR;

// 	previous time step part
	for(size_t i=0; i < m_prevSteps; ++i)
	{
		if(this->m_pDomDisc->assemble_defect
				(d, *(*m_pSolOld)[i], dofDistr, (*m_pTimeOld)[i],
				 s_m[i+1], s_a[i+1]*m_dt) != IAssemble_OK)
			return IAssemble_ERROR;
	}

//	we're done
	return IAssemble_OK;
}

template <typename TDoFDistribution, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_solution(vector_type& u, const dof_distribution_type& dofDistr)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:assemble_solution':"
				" Domain Discretization not set.\n");
		return IAssemble_ERROR;
	}

//	result
	IAssembleReturn res;

//	assemble solution
	res = this->m_pDomDisc->assemble_solution(u, dofDistr, m_futureTime);

//	interprete result
	switch(res)
	{
	case IAssemble_ERROR:
			UG_LOG("ERROR in assemble_solution.\n");
			return IAssemble_ERROR;
	case IAssemble_NOT_IMPLEMENTED:
			UG_LOG("ERROR in assemble_solution: function not implemented.\n");
			return IAssemble_ERROR;
	case IAssemble_TIME_INDEPENDENT:
			UG_LOG("ERROR in assemble_solution: Problem time independent.\n");
			return IAssemble_ERROR;

	case IAssemble_OK:
	case IAssemble_NON_LINEAR:	return IAssemble_OK;
	}

//	we're done
	return IAssemble_OK;
}

template <typename TDoFDistribution, typename TAlgebra >
IAssembleReturn
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, const vector_type& u,
                const dof_distribution_type& dofDistr)
{
	return IAssemble_NOT_IMPLEMENTED;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__THETA_TIME_STEP_IMPL__ */

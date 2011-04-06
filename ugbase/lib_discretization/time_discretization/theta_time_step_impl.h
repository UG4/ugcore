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
prepare_step(SolutionTimeSeries<vector_type>& prevSol,
             number dt)
{
//	perform checks
	if(prevSol.size() != m_prevSteps)
	{
		UG_LOG("ERROR in ThetaTimeDiscretization::prepare_step:"
			" Number of previous solutions must be "<<m_prevSteps<<".\n");
		return false;
	}
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:prepare_step': "
				"Domain Discretization not set.\n");
		return false;
	}

//	remember old values
	m_pPrevSol = &prevSol;

//	remember time step size
	m_dt = dt;

//	compute future time
	m_futureTime = m_dt + m_pPrevSol->time(0);

//	done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u,
                  const dof_distribution_type& dofDistr)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:assemble_jacobian':"
				" Domain Discretization not set.\n");
		return false;
	}

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	m_pPrevSol->push(*const_cast<vector_type*>(&u), m_futureTime);

//	assemble jacobian using current iterate
	if(this->m_pDomDisc->assemble_jacobian
			(J, u, m_futureTime, (*m_pPrevSol), dofDistr, s_m[0], s_a[0]*m_dt) != true)
		return false;

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();

//	we're done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u,
                const dof_distribution_type& dofDistr)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:assemble_defect':"
				" Domain Discretization not set.\n");
		return false;
	}


//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	m_pPrevSol->push(*const_cast<vector_type*>(&u), m_futureTime);

// 	future solution part
	if(this->m_pDomDisc->assemble_defect
			(d, u, m_futureTime, (*m_pPrevSol), dofDistr, s_m[0], s_a[0]*m_dt) != true)
		return false;

// 	previous time step part
	for(size_t i=1; i < m_pPrevSol->size(); ++i)
	{
		if(this->m_pDomDisc->assemble_defect
				(d,  m_pPrevSol->solution(i), m_pPrevSol->time(i), (*m_pPrevSol), dofDistr,
				 s_m[i], s_a[i]*m_dt) != true)
			return false;
	}

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();

//	we're done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_solution(vector_type& u, const dof_distribution_type& dofDistr)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:assemble_solution':"
				" Domain Discretization not set.\n");
		return false;
	}

//	result
	bool res;

//	assemble solution
	res = this->m_pDomDisc->assemble_solution(u, m_futureTime, dofDistr);

//	we're done
	return res;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, const vector_type& u,
                const dof_distribution_type& dofDistr)
{
	return false;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISCRETIZATION__TIME_DISCRETIZATION__THETA_TIME_STEP_IMPL__ */

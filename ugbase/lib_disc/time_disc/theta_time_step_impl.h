/*
 * theta_time_step_impl.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__
#define __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__

#include "theta_time_step.h"

namespace ug{

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
prepare_step(VectorTimeSeries<vector_type>& prevSol,
             number dt)
{
//	perform checks
	if(prevSol.size() < m_prevSteps)
	{
		UG_LOG("ERROR in ThetaTimeDiscretization::prepare_step:"
			" Number of previous solutions must at least "<<m_prevSteps<<".\n");
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

//	update scalings
	update_scaling(m_dt);

//	done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u,
                  const dof_distribution_type& dd)
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

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	J.resize(0,0);
	J.resize(numIndex, numIndex);
	J.set(0.0);

//	assemble jacobian using current iterate
	if(!this->m_pDomDisc->assemble_jacobian
			(J, *m_pPrevSol, m_vScaleStiff[0], dd))
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
                const dof_distribution_type& dd)
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

//	reset matrix to zero and resize
	const size_t numIndex = dd.num_indices();
	d.resize(numIndex);
	d.set(0.0);

	UG_ASSERT(m_pPrevSol->size() >= m_prevSteps + 1, "Wrong number of solutions")

// 	future solution part
	if(!this->m_pDomDisc->assemble_defect
			(d, *m_pPrevSol, m_vScaleMass, m_vScaleStiff, dd))
		return false;

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();

//	we're done
	return true;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, const dof_distribution_type& dd)
{
//	check domain disc
	if(this->m_pDomDisc == NULL)
	{
		UG_LOG("ERROR in 'ThetaTimeDiscretization:adjust_solution':"
				" Domain Discretization not set.\n");
		return false;
	}

//	result
	bool res;

//	assemble solution
	res = this->m_pDomDisc->adjust_solution(u, m_futureTime, dd);

//	we're done
	return res;
}

template <typename TDoFDistribution, typename TAlgebra >
bool
ThetaTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b,
                const dof_distribution_type& dd)
{
	return false;
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__ */

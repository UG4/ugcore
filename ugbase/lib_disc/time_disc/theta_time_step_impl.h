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
void
MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>::
prepare_step(VectorTimeSeries<vector_type>& prevSol,
             number dt)
{
//	perform checks
	if(prevSol.size() < m_prevSteps)
		UG_THROW_FATAL("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must at least "<<m_prevSteps<<".\n");

//	remember old values
	m_pPrevSol = &prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              *m_pPrevSol);
}

template <typename TDoFDistribution, typename TAlgebra >
void
MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u,
                  const dof_distribution_type& dd)
{
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
	try{
		this->m_rDomDisc.assemble_jacobian(J, *m_pPrevSol, m_vScaleStiff[0], dd);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TDoFDistribution, typename TAlgebra >
void
MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u,
                const dof_distribution_type& dd)
{
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
	try{
	this->m_rDomDisc.assemble_defect(d, *m_pPrevSol, m_vScaleMass, m_vScaleStiff, dd);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble defect.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TDoFDistribution, typename TAlgebra >
void
MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>::
adjust_solution(vector_type& u, const dof_distribution_type& dd)
{
//	assemble solution
	try{
	this->m_rDomDisc.adjust_solution(u, m_futureTime, dd);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot adjust solution.");
}

template <typename TDoFDistribution, typename TAlgebra >
void
MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b,
                const dof_distribution_type& dd)
{
	UG_THROW_FATAL("Not implemented.");
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__ */

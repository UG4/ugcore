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

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
prepare_step(VectorTimeSeries<vector_type>& prevSol,
             number dt)
{
//	perform checks
	if(prevSol.size() < m_prevSteps)
		UG_THROW_FATAL("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol.size() << " passed.\n");

//	remember old values
	m_pPrevSol = &prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              *m_pPrevSol);
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
prepare_step_elem(VectorTimeSeries<vector_type>& prevSol,
                  number dt, GridLevel gl)
{
//	perform checks
	if(prevSol.size() < m_prevSteps)
		UG_THROW_FATAL("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol.size() << " passed.\n");

//	remember old values
	m_pPrevSol = &prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              *m_pPrevSol);

// 	prepare timestep
	try{
		this->m_rDomDisc.prepare_timestep(*m_pPrevSol, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot prepare timestep.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u, GridLevel gl)
{
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW_FATAL("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.\n");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	m_pPrevSol->push(*const_cast<vector_type*>(&u), m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_rDomDisc.assemble_jacobian(J, *m_pPrevSol, m_vScaleStiff[0], gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u, GridLevel gl)
{
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW_FATAL("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.\n");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	m_pPrevSol->push(*const_cast<vector_type*>(&u), m_futureTime);

// 	future solution part
	try{
		this->m_rDomDisc.assemble_defect(d, *m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble defect.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
adjust_solution(vector_type& u, GridLevel gl)
{
//	assemble solution
	try{
		this->m_rDomDisc.adjust_solution(u, m_futureTime, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot adjust solution.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, GridLevel gl)
{
	UG_THROW_FATAL("Not implemented.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
finish_step_elem(VectorTimeSeries<vector_type>& prevSol,
                 number dt, GridLevel gl)
{
//	perform checks
	if(prevSol.size() < m_prevSteps)
		UG_THROW_FATAL("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol.size() << " passed.\n");

//	remember old values and values of current timestep
	m_pPrevSol = &prevSol;

// 	finish timestep
	try{
		this->m_rDomDisc.finish_timestep(*m_pPrevSol, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot finish timestep.");
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__ */

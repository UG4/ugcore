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
prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
             number dt)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_prepare_step, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol->size() << " passed.\n");

//	remember old values
	m_pPrevSol = prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              m_pPrevSol);
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
                  number dt, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_step_elem, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol->size() << " passed.\n");

//	remember old values
	m_pPrevSol = prevSol;

//	remember time step size
	m_dt = dt;

//	update scalings
	m_futureTime = update_scaling(m_vScaleMass, m_vScaleStiff,
	                              m_dt, m_pPrevSol->time(0),
	                              m_pPrevSol);

// 	prepare timestep
	try{
		this->m_spDomDisc->prepare_timestep(m_pPrevSol, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot prepare timestep.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_jacobian, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::assemble_jacobian:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
	m_pPrevSol->push(pU, m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_jacobian(J, m_pPrevSol, m_vScaleStiff[0], gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_defect, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::assemble_defect:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
	m_pPrevSol->push(pU, m_futureTime);

// 	future solution part
	try{
		this->m_spDomDisc->assemble_defect(d, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble defect.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
adjust_solution(vector_type& u, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_adjust_solution, "discretization MultiStepTimeDiscretization");
//	assemble solution
	try{
		this->m_spDomDisc->adjust_solution(u, m_futureTime, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot adjust solution.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_linear, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::assemble_linear:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&b), &DummyRefCount);
	m_pPrevSol->push(pU, m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_linear(A, b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_rhs(vector_type& b, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_rhs, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::assemble_linear:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&b), &DummyRefCount);
	m_pPrevSol->push(pU, m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_rhs(b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_rhs(vector_type& b, const vector_type& u, GridLevel gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_rhs, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::assemble_linear:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
	m_pPrevSol->push(pU, m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_rhs(b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
                 number dt, GridLevel gl)
{
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("ThetaTimeStep::prepare_step:"
						" Number of previous solutions must be at least "<<
						m_prevSteps <<", but only "<< prevSol->size() << " passed.\n");

//	remember old values and values of current timestep
	m_pPrevSol = prevSol;

// 	finish timestep
	try{
		this->m_spDomDisc->finish_timestep(m_pPrevSol, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot finish timestep.");
}

} // end namespace ug

#endif /* __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__ */

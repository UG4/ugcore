/*
 * Copyright (c) 2010-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__
#define __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP_IMPL__

#include "theta_time_step.h"
#include "common/math/misc/math_constants.h"


namespace ug{

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
             number dt)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_prepare_step, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::prepare_step:"
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

//	prepare time step (elemDisc-wise)
	try
	{
		this->m_spDomDisc->prepare_timestep(m_pPrevSol, m_futureTime);
	}
	UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot prepare time step.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
                  number dt, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_step_elem, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(prevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::prepare_step_elem:"
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

//	prepare time step (elemDisc-wise)
	try
	{
		this->m_spDomDisc->prepare_timestep(m_pPrevSol, m_futureTime, gl);
	}
	UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot prepare time step.");

// 	prepare timestep element-wise
	try{
		this->m_spDomDisc->prepare_timestep_elem(m_pPrevSol, gl);
	}UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot prepare timestep element-wise.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_jacobian, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::assemble_jacobian:"
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
	}UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_defect, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::assemble_defect:"
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
adjust_solution(vector_type& u, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_adjust_solution, "discretization MultiStepTimeDiscretization");
//	adjust solution
	try{
		this->m_spDomDisc->adjust_solution(u, m_futureTime, gl);
	}UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot adjust solution.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_linear, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::assemble_linear:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");


//	push unknown solution to solution time series (not used, but formally needed)
	m_pPrevSol->push(m_pPrevSol->latest(), m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_linear(A, b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot assemble jacobian.");

//	pop unknown solution from solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_rhs(vector_type& b, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_rhs, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::assemble_linear:"
				" Number of previous solutions must be at least "<<
				m_prevSteps <<", but only "<< m_pPrevSol->size() << " passed.");

//	push unknown solution to solution time series (not used, but formally needed)
	m_pPrevSol->push(m_pPrevSol->latest(), m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_rhs(b, m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl);
	}UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot assemble jacobian.");

//	pop unknown solution from solution time series
	m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_assemble_rhs, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::assemble_linear:"
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

template<typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
calc_error(const vector_type& u, error_vector_type* u_vtk)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_calc_error, "discretization MultiStepTimeDiscretization");
//	perform checks
	if(m_pPrevSol->size() < m_prevSteps)
		UG_THROW("MultiStepTimeDiscretization::calc_error:"
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

//	assemble error estimators using current iterate
	try
	{
		GridLevel gl = GridLevel();	// new TOP-SURFACE grid level
		this->m_spDomDisc->calc_error(m_pPrevSol, m_vScaleMass, m_vScaleStiff, gl, u_vtk);
	}
	UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot assemble error estimators.");

//	pop unknown solution to solution time series
	m_pPrevSol->remove_latest();
}


template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
finish_step(SmartPtr<VectorTimeSeries<vector_type> > currSol)
{
//	perform checks whether 'currSol' is a solutionTimeSeries only with the new values
	if (currSol->time(0) != m_futureTime)
		UG_THROW("MultiStepTimeDiscretization::finish_step:"
				" The solution of the SolutionTimeSeries used in this function"
				" does not coincide with the current solution! ");

// 	finish timestep using the current solution
	try
	{
		this->m_spDomDisc->finish_timestep(currSol);
	}
	UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot finish timestep.");
}

template <typename TAlgebra>
void MultiStepTimeDiscretization<TAlgebra>::
finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
                 const GridLevel& gl)
{
//	perform checks whether 'currSol' is a solutionTimeSeries only with the new values
	if(currSol->time(0) != m_futureTime)
		UG_THROW("MultiStepTimeDiscretization::finish_step_elem:"
				" The solution of the SolutionTimeSeries used in this function"
				" does not coincide with the current solution! ");

// 	finish timestep using the current solution
	try
	{
		this->m_spDomDisc->finish_timestep(currSol, gl);
	}
	UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot finish timestep.");

// 	finish timestep element-wise using the current solution
	try{
		this->m_spDomDisc->finish_timestep_elem(currSol, gl);
	}UG_CATCH_THROW("MultiStepTimeDiscretization: Cannot finish timestep element-wise.");
}

//////////////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////

template <typename TAlgebra>
void SDIRK<TAlgebra>::set_stage(size_t stage)
{
	m_stage = stage;
}

template <typename TAlgebra>
number SDIRK<TAlgebra>::update_scaling(std::vector<number>& vSM,
                                       std::vector<number>& vSA,
                                       number dt)
{
	if(m_order == 1) // Mittelpunkt
	{
		switch(m_stage)
		{
			case 1:
				vSM.resize(2);
				vSA.resize(2);
				vSM[0] = 1.;
				vSM[1] = -1.;
				vSA[0] = dt * 1./2.;
				vSA[1] = 0;
				return m_Time0 + 1./2. * dt;
			case 2:
				vSM.resize(3);
				vSA.resize(3);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = -1.;
				vSA[0] = 0;
				vSA[1] = dt;
				vSA[2] = 0;
				return m_Time0 + dt;
			default:
				UG_THROW("Midpoint scheme has only 2 stages")
		}
	}
	else if(m_order == 2) // Alexander, order 2
	{
		const number gamma = 1 - 1. / sqrt(2.);
		switch(m_stage)
		{
			case 1:
				vSM.resize(2);
				vSA.resize(2);
				vSM[0] = 1.;
				vSM[1] = -1.;
				vSA[0] = dt * gamma;
				vSA[1] = 0;
				return m_Time0 + gamma * dt;
			case 2:
				vSM.resize(3);
				vSA.resize(3);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = -1;
				vSA[0] = dt * gamma;
				vSA[1] = dt * (1. - gamma);
				vSA[2] = 0;
				return m_Time0 + dt;
			default:
				UG_THROW("Alexander(2) scheme has only 2 stages")
		}
	}
	else if(m_order == 3) // Alexander, order 3
	{
		const number alpha = 0.4358665215;// root of x^3-3x^2+3/2x-1/6 = 0
		const number tau = (1. + alpha)/2.;
		const number b1 = -(6*alpha*alpha-16*alpha+1.)/4.;
		const number b2 = (6*alpha*alpha-20*alpha+5.)/4.;
		switch(m_stage)
		{
			case 1:
				vSM.resize(2);
				vSA.resize(2);
				vSM[0] = 1.;
				vSM[1] = -1.;
				vSA[0] = dt * alpha;
				vSA[1] = 0;
				return m_Time0 + alpha * dt;
			case 2:
				vSM.resize(3);
				vSA.resize(3);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = -1;
				vSA[0] = dt * alpha;
				vSA[1] = dt * (tau-alpha);
				vSA[2] = 0;
				return m_Time0 + tau * dt;
			case 3:
				vSM.resize(4);
				vSA.resize(4);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = -1;
				vSA[0] = dt * alpha;
				vSA[1] = dt * b1;
				vSA[2] = dt * b2;
				vSA[3] = 0;
				return m_Time0 + dt;
			default:
				UG_THROW("Alexander(3) scheme has only 3 stages")
		}
	}
	else if(m_order == 4) // Hairer,Wanner, order 4, L-stable DIRK
	{
		switch(m_stage)
		{
			case 1:
				vSM.resize(2);
				vSA.resize(2);
				vSM[0] = 1.;
				vSM[1] = -1.;
				vSA[0] = dt * 1./4.;
				vSA[1] = 0;
				return m_Time0 + 1./4. * dt;
			case 2:
				vSM.resize(3);
				vSA.resize(3);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = -1;
				vSA[0] = dt * 1./4.;
				vSA[1] = dt * 1./2.;
				vSA[2] = 0;
				return m_Time0 + 3./4. * dt;
			case 3:
				vSM.resize(4);
				vSA.resize(4);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = -1;
				vSA[0] = dt * 1./4.;
				vSA[1] = dt * (-1./25.);
				vSA[2] = dt * 17./50.;
				vSA[3] = 0;
				return m_Time0 + 11./20. * dt;
			case 4:
				vSM.resize(5);
				vSA.resize(5);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = 0;
				vSM[4] = -1;
				vSA[0] = dt * 1./4.;
				vSA[1] = dt * (15./544.);
				vSA[2] = dt * (-137./2720.);
				vSA[3] = dt * (371./1360.);
				vSA[4] = 0;
				return m_Time0 + 1./2. * dt;
			case 5:
				vSM.resize(6);
				vSA.resize(6);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = 0;
				vSM[4] = 0;
				vSM[5] = -1;
				vSA[0] = dt * 1./4.;
				vSA[1] = dt * (-85./12.);
				vSA[2] = dt * (125./16.);
				vSA[3] = dt * (-49./48.);
				vSA[4] = dt * (25./24.);
				vSA[5] = 0;
				return m_Time0 + 1. * dt;
			default:
				UG_THROW("HairerWanner(4) scheme has only 5 stages")
		}
	}
	else if(m_order == 3) // Crouzeix, order 3
	{
		const number gamma = (3. + sqrt(3.))/6.;
		switch(m_stage)
		{
			case 1:
				vSM.resize(2);
				vSA.resize(2);
				vSM[0] = 1.;
				vSM[1] = -1.;
				vSA[0] = dt * gamma;
				vSA[1] = 0;
				return m_Time0 + gamma * dt;
			case 2:
				vSM.resize(3);
				vSA.resize(3);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = -1.;
				vSA[0] = dt * gamma;
				vSA[1] = dt * (1. - 2*gamma);
				vSA[2] = 0;
				return m_Time0 + (1-gamma)*dt;
			case 3:
				vSM.resize(4);
				vSA.resize(4);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = -1.;
				vSA[0] = 0;
				vSA[1] = dt * 1./2.;
				vSA[2] = dt * 1./2.;
				vSA[3] = 0;
				return m_Time0 + dt;
			default:
				UG_THROW("Crouzeix 3 scheme has only 3 stages")
		}
	}
	else if(m_order == 4) // Crouzeix, order 4
	{
		const number gamma = 1./2. + cos(M_PI/18.) / sqrt(3.);
		const number delta = 1./(6. * (2*gamma-1) * (2*gamma-1));
		switch(m_stage)
		{
			case 1:
				vSM.resize(2);
				vSA.resize(2);
				vSM[0] = 1.;
				vSM[1] = -1.;
				vSA[0] = dt * gamma;
				vSA[1] = 0;
				return m_Time0 + gamma * dt;
			case 2:
				vSM.resize(3);
				vSA.resize(3);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = -1.;
				vSA[0] = dt * gamma;
				vSA[1] = dt * (1./2. - gamma);
				vSA[2] = 0;
				return m_Time0 + 1./2.*dt;
			case 3:
				vSM.resize(4);
				vSA.resize(4);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = -1.;
				vSA[0] = dt * gamma;
				vSA[1] = dt * (1-4*gamma);
				vSA[2] = dt * 2*gamma;
				vSA[3] = 0;
				return m_Time0 + (1-gamma)* dt;
			case 4:
				vSM.resize(5);
				vSA.resize(5);
				vSM[0] = 1.;
				vSM[1] = 0;
				vSM[2] = 0;
				vSM[3] = 0;
				vSM[4] = -1.;
				vSA[0] = 0;
				vSA[1] = dt * delta;
				vSA[2] = dt * (1-2*delta);
				vSA[3] = dt * delta;
				vSA[4] = 0;
				return m_Time0 + dt;
			default:
				UG_THROW("Crouzeix 4 scheme has only 4 stages")
		}
	}
	else
		UG_THROW("SDIRK: "<< m_order <<" missing.");

}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
             number dt)
{
//	remember old values
	if(m_stage == 1){
		this->m_pPrevSol = SmartPtr<VectorTimeSeries<vector_type> >(
							new VectorTimeSeries<vector_type>);
		m_Time0 = prevSol->time(0);
		this->m_futureTime = m_Time0;
	}
	this->m_pPrevSol->push(prevSol->solution(0)->clone(), prevSol->time(0));

//	remember time step size
	this->m_dt = dt;

//	update scalings
	m_lastTime = this->m_futureTime;

	this->m_futureTime = update_scaling(this->m_vScaleMass, this->m_vScaleStiff,
	                                    this->m_dt);

//	prepare time step (elemDisc-wise)
	try
	{
		this->m_spDomDisc->prepare_timestep(this->m_pPrevSol, this->m_futureTime);
	}
	UG_CATCH_THROW("ThetaTimeStep: Cannot prepare time step.");
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl)
{
//	if(this->m_pPrevSol->size() < m_stage /*&& m_stage != num_stages()*/){
//		this->m_pPrevSol->push(u.clone(), m_lastTime);
//	}

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
	this->m_pPrevSol->push(pU, this->m_futureTime);

//	assemble jacobian using current iterate
	try{
		this->m_spDomDisc->assemble_jacobian(J, this->m_pPrevSol, this->m_vScaleStiff[0], gl);
	}UG_CATCH_THROW("SDIRK: Cannot assemble jacobian.");

//	pop unknown solution to solution time series
	this->m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl)
{
//	if(this->m_pPrevSol->size() < m_stage /*&& m_stage != num_stages()*/){
//		this->m_pPrevSol->push(u.clone(), m_lastTime);
//	}

//	push unknown solution to solution time series
//	ATTENTION: Here, we must cast away the constness of the solution, but note,
//			   that we pass pPrevSol as a const object in assemble_... Thus,
//			   the solution will not be changed there and we pop it from the
//			   Solution list afterwards, such that nothing happens to u
	// \todo: avoid this hack, use smart ptr properly
	int DummyRefCount = 2;
	SmartPtr<vector_type> pU(const_cast<vector_type*>(&u), &DummyRefCount);
	this->m_pPrevSol->push(pU, this->m_futureTime);

// 	future solution part
	try{
		this->m_spDomDisc->assemble_defect(d,  this->m_pPrevSol,  this->m_vScaleMass,  this->m_vScaleStiff, gl);
	}UG_CATCH_THROW("SDIRK: Cannot assemble defect.");

//	pop unknown solution to solution time series
	this->m_pPrevSol->remove_latest();
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
adjust_solution(vector_type& u, const GridLevel& gl)
{
	PROFILE_BEGIN_GROUP(MultiStepTimeDiscretization_adjust_solution, "discretization MultiStepTimeDiscretization");
//	adjust solution
	try{
		this->m_spDomDisc->adjust_solution(u, this->m_futureTime, gl);
	}UG_CATCH_THROW("ThetaTimeStep: Cannot adjust solution.");
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl)
{
	UG_THROW("Not implemented")
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
assemble_rhs(vector_type& b, const GridLevel& gl)
{
	UG_THROW("Not implemented")
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl)
{
	UG_THROW("Not implemented")
}

/* Please overwrite any of the following methods, if applicable:
template <typename TAlgebra>
void SDIRK<TAlgebra>::
prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
                  number dt, const GridLevel& gl)
{
	UG_THROW("Not implemented")
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
finish_step(SmartPtr<VectorTimeSeries<vector_type> > currSol)
{
	UG_THROW("Not implemented")
}

template <typename TAlgebra>
void SDIRK<TAlgebra>::
finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
                 const GridLevel& gl)
{
	UG_THROW("Not implemented")
}
*/


} // end namespace ug


#endif
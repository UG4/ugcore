/*
 * timestep.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__TIMESTEP__
#define __H__LIB_DISCRETIZATION__TIMESTEP__

#include "spacialdiscretization.h"
#include "numericalsolution.h"
#include "../common/types.h"
#include "newton.h"

namespace ug{

class TimeStep : public NonLinearAssemble {

	public:
		TimeStep(number theta) // theta = 0 -> Backward Euler
		{
			s_a = new number[2];
			s_m = new number[2];
			s_a[0] = 1.-theta;
			s_a[1] = theta;
			s_m[0] = 1.;
			s_m[1] = -1.;
			m_previousSteps = 1;
		}

		bool init(std::vector<NumericalSolution*>& u_old, std::vector<number>& time_old, number time_future, SpacialDiscretization& sd)
		{
			m_u_old = &u_old;
			m_time_old = &time_old;
			m_sd = &sd;
			m_time_future = time_future;

			return true;
		}

		bool assemble_defect(Vector& vec, NumericalSolution& u)
		{
			if(m_u_old->size() != m_previousSteps) return false;
			if(m_time_old->size() != m_previousSteps) return false;

			DirichletValues dirVal;

			number dt = m_time_future - (*m_time_old)[0];

			// future solution part
			if(m_sd->assemble_defect(vec, u, m_time_future, s_m[0], s_a[0]*dt) == false) return false;

			// previous time step part
			for(uint i=0; i < m_previousSteps; ++i)
			{
				if(m_sd->assemble_defect(vec, *(*m_u_old)[i], (*m_time_old)[i], s_m[i+1], s_a[i+1]*dt) == false) return false;
			}

			if(m_sd->get_dirichlet_values(u, dirVal) == false) return false;
			dirVal.set_zero_values(vec);

			return true;
		}

		bool assemble_jacobian(Matrix& mat, NumericalSolution& u)
		{
			if(m_u_old->size() != m_previousSteps) return false;
			if(m_time_old->size() != m_previousSteps) return false;

			DirichletValues dirVal;

			number dt = m_time_future - (*m_time_old)[0];

			// future solution part
			if(m_sd->assemble_jacobian(mat, u, m_time_future, s_m[0], s_a[0]*dt) == false) return false;

			if(m_sd->get_dirichlet_values(u, dirVal) == false) return false;
			dirVal.set_rows(mat);

			return true;
		}

		bool assemble_solution(NumericalSolution& u)
		{
			DirichletValues dirVal;

			if(m_sd->get_dirichlet_values(u, dirVal) == false) return false;

			dirVal.set_values(*u.GridVector());

			return true;
		}


		~TimeStep()
		{
			delete s_a;
			delete s_m;
		}

	private:
		int m_previousSteps;
		number *s_a;
		number *s_m;
		std::vector<NumericalSolution*> *m_u_old;
		std::vector<number> *m_time_old;
		SpacialDiscretization *m_sd;
		number m_time_future;

};




}



#endif /* __H__LIB_DISCRETIZATION__TIMESTEP__ */

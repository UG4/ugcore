/*
 * timestep.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__TIMESTEP__
#define __H__LIB_DISCRETIZATION__TIMESTEP__

// extern libraries
#include <deque>

// other ug libraries
#include "../common/types.h"

// modul intern libraries
#include "numericalsolution.h"
#include "assemble.h"

namespace ug{



/// Time Discretization Interface
/**
 * implements the time discretization.
 *
 * This class uses a ISpacialDiscratization in order to implement the IAssemble interface.
 *
 * After the method prepare step has been called, Jacobian/Defect can be computed.
 *
 */
template <int dim>
class ITimeDiscretization : public IAssemble<dim> {
	public:
		/// constructor
		/**
		 * \param[in] isd	SpacialDiscretization
		 */
		ITimeDiscretization(ISpacialDiscretization<dim>& isd)
		{
			_isd = &isd;
		}

		/// prepares the assembling of Defect/Jacobian (resp. Jacobian/rhs) for a time - dependent and nonlinear (resp. linear) problem
		/**
		 *	This function supplies the TimeDiscretization with previous time steps and step size before the assembling routines
		 *	can be called.
		 *
		 * \param[in] u_old 	a std::vector containing the solution at the previous time steps
		 * \param[in] time_old	a std::vector containing the time at the previous time steps
		 * \param[in] dt		size of time step
		 */
		virtual bool prepare_step(std::deque<NumericalSolution<dim>*>& u_old, std::deque<number>& time_old, number dt) = 0;

	protected:
		ISpacialDiscretization<dim>* _isd;

};


template <int dim>
class ThetaTimeDiscretization : public ITimeDiscretization<dim> {

	public:
		ThetaTimeDiscretization(ISpacialDiscretization<dim>& sd, number theta) : ITimeDiscretization<dim>(sd) // theta = 0 -> Backward Euler
		{
			s_a[0] = 1.-theta;
			s_a[1] = theta;
			s_m[0] = 1.;
			s_m[1] = -1.;
			_previousSteps = 1;
		}

		bool prepare_step(std::deque<NumericalSolution<dim>*>& u_old, std::deque<number>& time_old, number dt)
		{
			if(u_old.size() != _previousSteps)
			{
				std::cout << "ERROR in prepare_step: Number of previous solutions is not adequate for this time solver" << std::endl;
				return false;
			}
			if(time_old.size() != _previousSteps)
			{
				std::cout << "ERROR in prepare_step: Number of previous time steps is not adequate for this time solver" << std::endl;
				return false;
			}
			if(dt < 0.0)
			{
				std::cout << "ERROR in prepare_step: Time step size can not be negative." << std::endl;
				return false;
			}

			_u_old = &u_old;
			_time_old = &time_old;
			_dt = dt;
			_time_future = _dt + (*_time_old)[0];

			return true;
		}


		IAssembleReturn assemble_jacobian_defect(Matrix& J, Vector& d, const NumericalSolution<dim>& u, uint level = 0)
		{
			// future solution part
			if(this->_isd->assemble_defect(d, u, _time_future, s_m[0], s_a[0]*_dt, level) != IAssemble_OK) return IAssemble_ERROR;
			if(this->_isd->assemble_jacobian(J, u, _time_future, s_m[0], s_a[0]*_dt, level) != IAssemble_OK) return IAssemble_ERROR;

			// previous time step part
			for(uint i=0; i < _previousSteps; ++i)
			{
				if(this->_isd->assemble_defect(d, *(*_u_old)[i], (*_time_old)[i], s_m[i+1], s_a[i+1]*_dt, level) == IAssemble_OK) return IAssemble_ERROR;
			}

			return IAssemble_OK;
		}

		IAssembleReturn assemble_jacobian(Matrix& J, NumericalSolution<dim>& u, uint level = 0)
		{
			if(this->_isd->assemble_jacobian(J, u, _time_future, s_m[0], s_a[0]*_dt, level) != IAssemble_OK) return IAssemble_ERROR;

			return IAssemble_OK;
		}

		IAssembleReturn assemble_defect(Vector& d, NumericalSolution<dim>& u, uint level = 0)
		{
			// future solution part
			if(this->_isd->assemble_defect(d, u, _time_future, s_m[0], s_a[0]*_dt, level) != IAssemble_OK) return IAssemble_ERROR;

			// previous time step part
			for(uint i=0; i < _previousSteps; ++i)
			{
				if(this->_isd->assemble_defect(d, *(*_u_old)[i], (*_time_old)[i], s_m[i+1], s_a[i+1]*_dt, level) != IAssemble_OK) return IAssemble_ERROR;
			}

			return IAssemble_OK;
		}

		IAssembleReturn assemble_solution(NumericalSolution<dim>& u, uint level = 0)
		{
			IAssembleReturn res;

			res = this->_isd->assemble_solution(u, _time_future, level);

			switch(res)
			{
			case IAssemble_ERROR:
					std::cout << "ERROR in assemble_solution" << std::endl; return IAssemble_ERROR;
			case IAssemble_NOT_IMPLEMENTED:
					std::cout << "ERROR in assemble_solution: function not implemented" << std::endl; return IAssemble_ERROR;
			case IAssemble_TIME_INDEPENDENT:
					std::cout << "ERROR in assemble_solution: Problem time independent" << std::endl; return IAssemble_ERROR;

			case IAssemble_OK:
			case IAssemble_NON_LINEAR:	return IAssemble_OK;
			}
			return IAssemble_OK;
		}

		IAssembleReturn assemble_linear(Matrix& A, Vector& b, uint level = 0)
		{
			return IAssemble_NOT_IMPLEMENTED;
		}

	private:
		uint _previousSteps;
		number s_a[2];
		number s_m[2];
		std::deque<NumericalSolution<dim>*> *_u_old;
		std::deque<number> *_time_old;
		number _dt; // current time step
		number _time_future;

};




}



#endif /* __H__LIB_DISCRETIZATION__TIMESTEP__ */

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
#include "common/common.h"

// modul intern libraries
#include "lib_discretization/assemble.h"

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
template <typename TAlgebra, typename TDiscreteFunction>
class ITimeDiscretization : public IAssemble<TAlgebra, TDiscreteFunction> {
	public:
		typedef TDiscreteFunction discrete_function_type;
		typedef TAlgebra algebra_type;

	public:
		/// constructor
		/**
		 * \param[in] isd	SpacialDiscretization
		 */
		ITimeDiscretization(IDomainDiscretization<algebra_type, discrete_function_type>& dd) : m_dd(dd)
		{}

		/// prepares the assembling of Defect/Jacobian (resp. Jacobian/rhs) for a time - dependent and nonlinear (resp. linear) problem
		/**
		 *	This function supplies the TimeDiscretization with previous time steps and step size before the assembling routines
		 *	can be called.
		 *
		 * \param[in] u_old 	a std::vector containing the solution at the previous time steps
		 * \param[in] time_old	a std::vector containing the time at the previous time steps
		 * \param[in] dt		size of time step
		 */
		virtual bool prepare_step(std::deque<discrete_function_type*>& u_old, std::deque<number>& time_old, number dt) = 0;

		std::size_t num_fct() const
		{
			return m_dd.num_fct();
		}

		bool boundary_value(number& val, typename discrete_function_type::position_type& corner, uint fct, number time = 0.0)
		{
			return m_dd.boundary_value(val, corner, fct, time);
		}

	protected:
		IDomainDiscretization<algebra_type, discrete_function_type>& m_dd;

};


template <typename TAlgebra, typename TDiscreteFunction>
class ThetaTimeDiscretization : public ITimeDiscretization<TAlgebra, TDiscreteFunction> {
	public:
		typedef TDiscreteFunction discrete_function_type;

		// type of algebra
		typedef TAlgebra algebra_type;

		// type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

		// type of algebra vector
		typedef typename algebra_type::vector_type vector_type;


	public:
		ThetaTimeDiscretization(IDomainDiscretization<algebra_type, discrete_function_type>& sd, number theta) : ITimeDiscretization<TAlgebra, TDiscreteFunction>(sd) // theta = 0 -> Backward Euler
		{
			s_a[0] = 1.-theta;
			s_a[1] = theta;
			s_m[0] = 1.;
			s_m[1] = -1.;
			m_previousSteps = 1;
		}

		bool prepare_step(std::deque<discrete_function_type*>& u_old, std::deque<number>& time_old, number dt)
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


		IAssembleReturn assemble_jacobian_defect(matrix_type& J, vector_type& d, const discrete_function_type& u)
		{
			// future solution part
			if(this->m_dd.assemble_defect(d, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;
			if(this->m_dd.assemble_jacobian(J, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;

			// previous time step part
			for(uint i=0; i < m_previousSteps; ++i)
			{
				if(this->m_dd.assemble_defect(d, *(*m_u_old)[i], (*m_time_old)[i], s_m[i+1], s_a[i+1]*m_dt) == IAssemble_OK) return IAssemble_ERROR;
			}

			return IAssemble_OK;
		}

		IAssembleReturn assemble_jacobian(matrix_type& J, const discrete_function_type& u)
		{
			if(this->m_dd.assemble_jacobian(J, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;

			return IAssemble_OK;
		}

		IAssembleReturn assemble_defect(vector_type& d, const discrete_function_type& u)
		{
			// future solution part
			if(this->m_dd.assemble_defect(d, u, m_time_future, s_m[0], s_a[0]*m_dt) != IAssemble_OK) return IAssemble_ERROR;

			// previous time step part
			for(uint i=0; i < m_previousSteps; ++i)
			{
				if(this->m_dd.assemble_defect(d, *(*m_u_old)[i], (*m_time_old)[i], s_m[i+1], s_a[i+1]*m_dt) != IAssemble_OK) return IAssemble_ERROR;
			}

			return IAssemble_OK;
		}

		IAssembleReturn assemble_solution(discrete_function_type& u)
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

		IAssembleReturn assemble_linear(matrix_type& A, vector_type& b, const discrete_function_type& u)
		{
			return IAssemble_NOT_IMPLEMENTED;
		}

	private:
		uint m_previousSteps;
		number s_a[2];
		number s_m[2];
		std::deque<discrete_function_type*> *m_u_old;
		std::deque<number> *m_time_old;
		number m_dt; // current time step
		number m_time_future;

};




}



#endif /* __H__LIB_DISCRETIZATION__TIMESTEP__ */

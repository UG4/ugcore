/*
 * time_discretization_interface.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__TIME_DISCRETIZATION_INTERFACE__
#define __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__TIME_DISCRETIZATION_INTERFACE__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// modul intern libraries
#include "lib_discretization/spacial_discretization/domain_discretization_interface.h"

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
template <	typename TDiscreteFunction,
			typename TAlgebra = typename TDiscreteFunction::algebra_type >
class ITimeDiscretization : public IAssemble<TDiscreteFunction, TAlgebra> {
	public:
		typedef TDiscreteFunction discrete_function_type;
		typedef TAlgebra algebra_type;

	public:
		/// constructor
		/**
		 * \param[in] isd	SpacialDiscretization
		 */
		ITimeDiscretization(IDomainDiscretization<discrete_function_type, algebra_type>& dd) : m_dd(dd)
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

		// returns number of previous time steps needed (i.e. size of deque for time_old and u_old)
		virtual size_t num_prev_steps() = 0;

		size_t num_fct() const{return m_dd.num_fct();}

		virtual bool is_dirichlet(int si, size_t fct) {return m_dd.is_dirichlet(si, fct);}

	protected:
		IDomainDiscretization<discrete_function_type, algebra_type>& m_dd;

};

} // end namespace ug

#endif /* __H__LIB_DISCRETIZATION__TIME_DISRETIZATION__TIME_DISCRETIZATION_INTERFACE__ */

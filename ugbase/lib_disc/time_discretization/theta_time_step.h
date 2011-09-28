/*
 * theta_time_step.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__TIME_DISCRETIZATION__THETA_TIME_STEP__
#define __H__UG__LIB_DISC__TIME_DISCRETIZATION__THETA_TIME_STEP__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// modul intern libraries
#include "lib_disc/time_discretization/time_discretization_interface.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// theta time stepping scheme
/**
 * This time stepping scheme discretizes equations of the form
 * \f[
 * 	\partial_t u(t) = f(t)
 * \f]
 * as
 * \f[
 * 	\frac{u(t^{k+1}) - u(t^k)}{\Delta t} = (1-\theta) \cdot f(t^{k+1})
 * 												+ \theta \cdot f(t^k)
 * \f]
 *
 * Thus, for \f$\theta = 1 \f$ this is the Backward-Euler time stepping.
 */
template <	typename TDoFDistribution,
			typename TAlgebra>
class ThetaTimeDiscretization
	: public ITimeDiscretization<TDoFDistribution, TAlgebra>
{
	public:
	//	DoF Distribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	// 	Type of algebra
		typedef TAlgebra algebra_type;

	// 	Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	// 	Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	// 	Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, algebra_type>
			domain_discretization_type;

	public:
		// theta = 0 -> Backward Euler
		ThetaTimeDiscretization(domain_discretization_type& sd, number theta)
			: ITimeDiscretization<TDoFDistribution, TAlgebra>(sd),
			  m_pPrevSol(NULL)
		{
			set_theta(1.0);
		}

		ThetaTimeDiscretization(){set_theta(1.0);}

	///	sets the theta value
		void set_theta(number theta)
		{
			m_theta = theta;
			m_prevSteps = 1;
		}

	///	sets the domain discretization
		void set_domain_discretization(domain_discretization_type& dd)
		{
			this->m_pDomDisc = &dd;
		}

	/// \copydoc ITimeDiscretization::num_prev_steps()
		size_t num_prev_steps() {return m_prevSteps;}

	///	\copydoc ITimeDiscretization::prepare_step()
		virtual bool prepare_step(VectorTimeSeries<vector_type>& prevSol,
		                          number dt);

	public:
	//	Implements the assemble interface
		bool assemble_jacobian(matrix_type& J,
		                       const vector_type& u,
		                       const dof_distribution_type& dofDistr);

		bool assemble_defect(vector_type& d,
		                     const vector_type& u,
		                     const dof_distribution_type& dofDistr);

		bool assemble_solution(vector_type& u,
		                       const dof_distribution_type& dofDistr);

		bool assemble_linear(matrix_type& A, vector_type& b,
		                     const vector_type& u,
		                     const dof_distribution_type& dofDistr);

	private:
	///	updates the scaling factors
		void update_scaling(number dt)
		{
			m_vScaleMass.resize(2);
			m_vScaleMass[0] = 1.;
			m_vScaleMass[1] = -1.;

			m_vScaleStiff.resize(2);
			m_vScaleStiff[0] = (1.- m_theta) * dt;
			m_vScaleStiff[1] = m_theta * dt;
		}


		number m_theta;

		size_t m_prevSteps;
		std::vector<number> m_vScaleMass;
		std::vector<number> m_vScaleStiff;

		VectorTimeSeries<vector_type>* m_pPrevSol;
		number m_dt; 							///< Time Step size
		number m_futureTime;					///< Future Time
};

} // end namespace ug

/// }

// include implementation
#include "theta_time_step_impl.h"

#endif /* __H__UG__LIB_DISC__TIME_DISCRETIZATION__THETA_TIME_STEP__ */

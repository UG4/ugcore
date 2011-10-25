/*
 * theta_time_step.h
 *
 *  Created on: 23.10.2009
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP__
#define __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP__

// extern libraries
#include <deque>

// other ug libraries
#include "common/common.h"

// modul intern libraries
#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/local_finite_element/common/lagrange1d.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// multi step time stepping scheme
template <	typename TDoFDistribution,
			typename TAlgebra>
class MultiStepTimeDiscretization
	: public ITimeDiscretization<TDoFDistribution, TAlgebra>
{
	public:
	///	DoF Distribution Type
		typedef IDoFDistribution<TDoFDistribution> dof_distribution_type;

	/// Type of algebra
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	/// Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, algebra_type>
			domain_discretization_type;

	public:
	/// constructor
		MultiStepTimeDiscretization(domain_discretization_type& sd)
			: ITimeDiscretization<TDoFDistribution, TAlgebra>(sd),
			  m_pPrevSol(NULL)
		{}

	/// \copydoc ITimeDiscretization::num_prev_steps()
		size_t num_prev_steps() {return m_prevSteps;}

	///	\copydoc ITimeDiscretization::prepare_step()
		virtual void prepare_step(VectorTimeSeries<vector_type>& prevSol,
		                          number dt);

	public:
		void assemble_jacobian(matrix_type& J, const vector_type& u,
		                       const dof_distribution_type& dd);

		void assemble_defect(vector_type& d, const vector_type& u,
		                     const dof_distribution_type& dd);

		void assemble_linear(matrix_type& A, vector_type& b,
		                     const dof_distribution_type& dd);

		void adjust_solution(vector_type& u,
		                       const dof_distribution_type& dd);

	protected:
	///	updates the scaling factors
		virtual void update_scaling(number dt) = 0;

		size_t m_prevSteps;					///< number of previous steps needed.
		std::vector<number> m_vScaleMass;	///< Scaling for mass part
		std::vector<number> m_vScaleStiff;	///< Scaling for stiffness part

		VectorTimeSeries<vector_type>* m_pPrevSol;	///< Previous solutions
		number m_dt; 								///< Time Step size
		number m_futureTime;						///< Future Time
};

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
	: public MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>
{
	public:
	///	Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, TAlgebra>
			domain_discretization_type;

	public:
	/// constructor
		ThetaTimeDiscretization(domain_discretization_type& sd)
			: MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>(sd)
		{
			set_theta(1.0);
			this->m_prevSteps = 1;
		}

	/// theta = 0 -> Backward Euler
		ThetaTimeDiscretization(domain_discretization_type& sd, number theta)
			: MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>(sd)
		{
			set_theta(theta);
			this->m_prevSteps = 1;
		}

	///	sets the theta value
		void set_theta(number theta) {m_theta = theta;}

	protected:
		virtual void update_scaling(number dt)
		{
			this->m_vScaleMass.resize(2);
			this->m_vScaleMass[0] = 1.;
			this->m_vScaleMass[1] = -1.;

			this->m_vScaleStiff.resize(2);
			this->m_vScaleStiff[0] = (1.- m_theta) * dt;
			this->m_vScaleStiff[1] = m_theta * dt;
		}

		number m_theta;
};


template <	typename TDoFDistribution,
			typename TAlgebra>
class BDF
	: public MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>
{
	public:
	///	Domain Discretization type
		typedef IDomainDiscretization<TDoFDistribution, TAlgebra>
			domain_discretization_type;

	public:
	/// constructor
		BDF(domain_discretization_type& sd)
			: MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>(sd)
		{
			set_order(2);
		}

	/// theta = 0 -> Backward Euler
		BDF(domain_discretization_type& sd, size_t order)
			: MultiStepTimeDiscretization<TDoFDistribution, TAlgebra>(sd)
		{
			set_order(order);
		}

	///	sets the theta value
		void set_order(size_t order) {m_order = order; this->m_prevSteps = order;}

	protected:
		virtual void update_scaling(number dt)
		{
		//	get reference to scaling factors
			std::vector<number>& vSM = this->m_vScaleMass;
			std::vector<number>& vSA = this->m_vScaleStiff;

		//	resize scaling factors
			vSM.resize(m_order+1);
			vSA.resize(m_order+1);

		//	get time points
			VectorTimeSeries<typename TAlgebra::vector_type>& prevSol= *this->m_pPrevSol;
			if(prevSol.size() < m_order+1)
				UG_THROW_FATAL("BDF("<<m_order<<") needs at least "<< m_order <<
				               " previous solutions, but only "<<prevSol.size()-1<<"passed.");

			std::vector<number> vTimePoint(m_order+1);
			for(size_t i = 0; i <= m_order; ++i)
				vTimePoint[i] = prevSol.time(i);

		//	create Lagrange Polynoms with given time steps
			for(size_t i = 0; i <= m_order; ++i)
			{
			//	create polynom
				Lagrange1D polynom(i, vTimePoint);

			//	get derivative
				Polynomial1D deriv = polynom.derivative();

			//	evaluate at future time point
				vSM[i] = deriv.value(vTimePoint[0]);
			}

		//	only first value of vSA != 0
			const number scale = 1.0 / vSM[0];
			vSA[0] = scale;
			for(size_t i = 1; i <= m_order; ++i) vSA[i] = 0;

		//	scale first s_m to 1.0
			for(size_t i = 0; i <= m_order; ++i) vSM[i] *= scale;
		}

		number m_order;
};


} // end namespace ug

/// }

// include implementation
#include "theta_time_step_impl.h"

#endif /* __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP__ */

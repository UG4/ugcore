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
#include <cmath>

// other ug libraries
#include "common/common.h"

// modul intern libraries
#include "lib_disc/time_disc/time_disc_interface.h"
#include "lib_disc/local_finite_element/common/lagrange1d.h"

namespace ug{

/// \ingroup lib_disc_time_assemble
/// @{

/// multi step time stepping scheme
template <typename TAlgebra>
class MultiStepTimeDiscretization
	: public ITimeDiscretization<TAlgebra>
{
	public:
	/// Type of algebra
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	/// Domain Discretization type
		typedef IDomainDiscretization<algebra_type>	domain_discretization_type;

	public:
	/// constructor
		MultiStepTimeDiscretization(SmartPtr<IDomainDiscretization<algebra_type> > spDD)
			: ITimeDiscretization<TAlgebra>(spDD),
			  m_pPrevSol(NULL)
		{}

		virtual ~MultiStepTimeDiscretization(){};

	/// \copydoc ITimeDiscretization::num_prev_steps()
		virtual size_t num_prev_steps() const {return m_prevSteps;}

	///	\copydoc ITimeDiscretization::prepare_step()
		virtual void prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
		                          number dt);

	///	\copydoc ITimeDiscretization::prepare_step_elem()
		virtual void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
		                               number dt, const GridLevel& gl);

	///	\copydoc ITimeDiscretization::finish_step_elem()
		virtual void finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
		                              const GridLevel& gl);

		virtual number future_time() const {return m_futureTime;}

	public:
		void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl);

		void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl);

		void assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl);

		void assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl);

		void assemble_rhs(vector_type& b, const GridLevel& gl);

		void adjust_solution(vector_type& u, const GridLevel& gl);

	protected:
	///	updates the scaling factors, returns the future time
		virtual number update_scaling(std::vector<number>& vSM,
		                              std::vector<number>& vSA,
		                              number dt, number currentTime,
		                              ConstSmartPtr<VectorTimeSeries<vector_type> > prevSol) = 0;

		size_t m_prevSteps;					///< number of previous steps needed.
		std::vector<number> m_vScaleMass;	///< Scaling for mass part
		std::vector<number> m_vScaleStiff;	///< Scaling for stiffness part

		SmartPtr<VectorTimeSeries<vector_type> > m_pPrevSol;	///< Previous solutions
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
template <typename TAlgebra>
class ThetaTimeStep
	: public MultiStepTimeDiscretization<TAlgebra>
{
	public:
	///	Domain Discretization type
		typedef IDomainDiscretization<TAlgebra> domain_discretization_type;

	/// Type of algebra
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// default constructor (implicit Euler)
		ThetaTimeStep(SmartPtr<IDomainDiscretization<TAlgebra> > spDD)
			: MultiStepTimeDiscretization<TAlgebra>(spDD),
			  m_stage(1), m_scheme("Theta")
		{
			set_theta(1.0);
			this->m_prevSteps = 1;
		}

	/// theta = 1.0 -> Implicit Euler, 0.0 -> Explicit Euler
		ThetaTimeStep(SmartPtr<IDomainDiscretization<TAlgebra> > spDD, number theta)
			: MultiStepTimeDiscretization<TAlgebra>(spDD),
			  m_stage(1), m_scheme("Theta")
		{
			set_theta(theta);
			this->m_prevSteps = 1;
		}

	/// theta = 1.0 -> Implicit Euler, 0.0 -> Explicit Euler
		ThetaTimeStep(SmartPtr<IDomainDiscretization<TAlgebra> > spDD, const char* scheme)
			: MultiStepTimeDiscretization<TAlgebra>(spDD),
			  m_stage(1), m_scheme(scheme)
		{
			set_theta(1.0);
			this->m_prevSteps = 1;
		}

		virtual ~ThetaTimeStep() {};

	///	sets the scheme
		void set_scheme(const char* scheme) {m_scheme = scheme;}

	///	returns number of stages
		virtual size_t num_stages() const
		{
			if		(m_scheme == "Theta") 		return 1;
			else if (m_scheme == "Alexander")	return 2;
			else if	(m_scheme == "FracStep") 	return 3;
			else UG_THROW("Step Scheme not recognized.");
		}

	///	sets the stage
		virtual void set_stage(size_t stage) {m_stage = stage;}

	///	sets the theta value
		void set_theta(number theta) {m_theta = theta;}

	protected:
		virtual number update_scaling(std::vector<number>& vSM,
		                              std::vector<number>& vSA,
		                              number dt, number currentTime,
		                              ConstSmartPtr<VectorTimeSeries<vector_type> > prevSol)
		{
		//	resize scaling factors
			vSM.resize(2);
			vSM[0] = 1.;
			vSM[1] = -1.;

			if(m_scheme == "Theta")
			{
				vSA.resize(2);
				vSA[0] = (m_theta) * dt;
				vSA[1] = (1.- m_theta) * dt;
				return currentTime + dt;
			}
			else if(m_scheme == "Alexander")
			{
				vSA.resize(2);
				const number gamma = 1 - 1. / sqrt(2.);
				switch(m_stage)
				{
					case 1:
						vSA[0] = gamma * dt;
						vSA[1] = 0;
						return currentTime + gamma * dt;
					case 2:
						vSA[0] = gamma * dt;
						vSA[1] = (1.- 2*gamma) * dt;
						return currentTime + (1 - gamma) * dt;
					default:
						UG_THROW("Alexander scheme has only 2 stages")
				}
			}
			else if(m_scheme == "FracStep")
			{
				vSA.resize(2);
				switch(m_stage)
				{
					case 1:
						vSA[0] = (2.-sqrt(2.)) 		* (1-1./sqrt(2.)) * dt;
						vSA[1] = (1.- (2.-sqrt(2.)))* (1-1./sqrt(2.)) * dt;
						return currentTime + (1-1./sqrt(2.)) * dt;
					case 2:
						vSA[0] = (sqrt(2.)-1) 		* (sqrt(2.)-1) * dt;
						vSA[1] = (1.- (sqrt(2.)-1)) * (sqrt(2.)-1) * dt;
						return currentTime + (sqrt(2.)-1) * dt;
					case 3:
						vSA[0] = (2.-sqrt(2.)) 		* (1-1./sqrt(2.)) * dt;
						vSA[1] = (1.- (2.-sqrt(2.)))* (1-1./sqrt(2.)) * dt;
						return currentTime + (1-1./sqrt(2.)) * dt;
					default:
						UG_THROW("Alexander scheme has only 2 stages")
				}
			}
			else
				UG_THROW("Unknown Multi-Stage Theta Scheme: "<< m_scheme<<".");

		}

		number m_theta;

		size_t m_stage;

		std::string m_scheme;
};


template <typename TAlgebra>
class BDF
	: public MultiStepTimeDiscretization<TAlgebra>
{
	public:
	///	Domain Discretization type
		typedef IDomainDiscretization<TAlgebra>	domain_discretization_type;

	/// Type of algebra
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// constructor
		BDF(SmartPtr<IDomainDiscretization<TAlgebra> > spDD)
			: MultiStepTimeDiscretization<TAlgebra>(spDD)
		{
			set_order(2);
		}

	/// theta = 0 -> Backward Euler
		BDF(SmartPtr<IDomainDiscretization<TAlgebra> > spDD, size_t order)
			: MultiStepTimeDiscretization<TAlgebra>(spDD)
		{
			set_order(order);
		}

		virtual ~BDF() {};

	///	sets the theta value
		void set_order(size_t order) {m_order = order; this->m_prevSteps = order;}

	///	returns the number of stages
		virtual size_t num_stages() const {return 1;}

	///	sets the stage
		virtual void set_stage(size_t stage)
		{
			if(stage!=1) UG_THROW("BDF has only one stage.");
		}

	protected:
		virtual number update_scaling(std::vector<number>& vSM,
		                              std::vector<number>& vSA,
		                              number dt, number currentTime,
		                              ConstSmartPtr<VectorTimeSeries<vector_type> > prevSol)
		{
		//	resize scaling factors
			vSM.resize(m_order+1);
			vSA.resize(m_order+1);

		//	future time
			const number futureTime = currentTime + dt;

		//	get time points
			if(prevSol->size() < m_order)
				UG_THROW("BDF("<<m_order<<") needs at least "<< m_order <<
				               " previous solutions, but only "<<prevSol->size()<<"passed.");

			std::vector<number> vTimePoint(m_order+1);
			vTimePoint[0] = futureTime;
			for(size_t i = 1; i <= m_order; ++i)
				vTimePoint[i] = prevSol->time(i-1);

		//	evaluate derivative of Lagrange Polynoms at future time
			for(size_t i = 0; i <= m_order; ++i)
			{
				vSM[i] = 0;

				for(size_t j = 0; j < vTimePoint.size(); ++j)
				{
					if(j == i) continue;

					number prod = 1;
					for(size_t k = 0; k < vTimePoint.size(); ++k)
					{
						if(k == i) continue;
						if(k == j) continue;

						prod *= (vTimePoint[0]-vTimePoint[k])/
								(vTimePoint[i]-vTimePoint[k]);
					}
					prod *= 1./(vTimePoint[i]-vTimePoint[j]);

					vSM[i] += prod;
				}
			}

		//	only first value of vSA != 0
			const number scale = 1.0 / vSM[0];
			vSA[0] = scale;
			for(size_t i = 1; i <= m_order; ++i) vSA[i] = 0;

		//	scale first s_m to 1.0
			for(int i = m_order; i >= 0; --i) vSM[i] *= scale;

			return futureTime;
		}

		size_t m_order;
};


/// Singly Diagonal Implicit Runge Kutta Method
template <typename TAlgebra>
class SDIRK
	: public MultiStepTimeDiscretization<TAlgebra>
{
	public:
	///	Domain Discretization type
		typedef IDomainDiscretization<TAlgebra> domain_discretization_type;

	/// Type of algebra
		typedef TAlgebra algebra_type;

	/// Type of algebra matrix
		typedef typename algebra_type::matrix_type matrix_type;

	/// Type of algebra vector
		typedef typename algebra_type::vector_type vector_type;

	public:
	/// default constructor (implicit Euler)
		SDIRK(SmartPtr<IDomainDiscretization<TAlgebra> > spDD)
			: MultiStepTimeDiscretization<TAlgebra>(spDD),
			  m_stage(1), m_order(1)
		{
			this->m_prevSteps = 1;
		}

	/// theta = 1.0 -> Implicit Euler, 0.0 -> Explicit Euler
		SDIRK(SmartPtr<IDomainDiscretization<TAlgebra> > spDD, int order)
			: MultiStepTimeDiscretization<TAlgebra>(spDD),
			  m_order(order)
		{
			set_order(m_order);
			this->m_prevSteps = 1;
		}

		virtual ~SDIRK() {};

	///	sets the scheme
		void set_order(int order) {
			if(order > 4) UG_THROW("DIRK: Only up to order 4.");
			m_order = order;
		}

	///	returns number of stages
		virtual size_t num_stages() const {
			switch(m_order){
				case 1: return 2;
				case 2: return 3;
				case 3: return 3;
				case 4: return 4;
				default: UG_THROW("DIRK: Only up to order 4.")
			}
		}

	///	sets the stage
		virtual void set_stage(size_t stage);

	public:
		virtual void prepare_step(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
								  number dt);
		virtual void prepare_step_elem(SmartPtr<VectorTimeSeries<vector_type> > prevSol,
									   number dt, const GridLevel& gl);
		virtual void finish_step_elem(SmartPtr<VectorTimeSeries<vector_type> > currSol,
									  const GridLevel& gl);

	public:
		void assemble_jacobian(matrix_type& J, const vector_type& u, const GridLevel& gl);

		void assemble_defect(vector_type& d, const vector_type& u, const GridLevel& gl);

		void assemble_linear(matrix_type& A, vector_type& b, const GridLevel& gl);

		void assemble_rhs(vector_type& b, const vector_type& u, const GridLevel& gl);

		void assemble_rhs(vector_type& b, const GridLevel& gl);

		void adjust_solution(vector_type& u, const GridLevel& gl);

	protected:
		virtual number update_scaling(std::vector<number>& vSM,
									  std::vector<number>& vSA,
									  number dt);

		virtual number update_scaling(std::vector<number>& vSM,
		                              std::vector<number>& vSA,
		                              number dt, number currentTime,
		                              ConstSmartPtr<VectorTimeSeries<vector_type> > prevSol) {
			UG_THROW("Not used.");
		}

		size_t m_stage;
		int m_order;
		number m_Time0;
		number m_lastTime;
};


} // end namespace ug


/// }

// include implementation
#include "theta_time_step_impl.h"

#endif /* __H__UG__LIB_DISC__TIME_DISC__THETA_TIME_STEP__ */

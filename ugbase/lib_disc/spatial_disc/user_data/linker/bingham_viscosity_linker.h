/*
 * Copyright (c) 2011-2015:  G-CSC, Goethe University Frankfurt
 * Author: Jonas Simon
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

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__BINGHAM_VISCOSITY_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__BINGHAM_VISCOSITY_LINKER__

#include "linker.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Bingham Viscosity linker
////////////////////////////////////////////////////////////////////////////////

/// Linker for the Bingham viscosity
/**
 * This linker computes the Bingham viscosity \f$ \mathbf{q} = - \frac{\mathbf{K}}{\mu} ( \nabla p - \rho \mathbf{g} ) \f$,
 * where
 * <ul>
 * <li> \f$ \mathbf{K} \f$		permeability
 * <li> \f$ \mu \f$				viscosity
 * <li> \f$ \rho \f$			density
 * <li> \f$ \mathbf{g} \f$		gravity
 * <li> \f$ \nabla p \f$		pressure gradient
 * </ul>
 * are input parameters.
 */

template <int dim>
class BinghamViscosityLinker
	: public StdDataLinker< BinghamViscosityLinker<dim>, number, dim>
{
	//	Base class type
	typedef StdDataLinker< BinghamViscosityLinker<dim>, number, dim> base_type;

	//  Constructor
	public:
		BinghamViscosityLinker() :
			m_spViscosity(NULL), m_spDViscosity(NULL),
			m_spDensity(NULL), m_spDDensity(NULL),
			m_spYieldStress(NULL), m_spDYieldStress(NULL),
			m_spVelocityGrad(NULL), m_spDVelocityGrad(NULL)
		{
		//	this linker needs exactly four input
			this->set_num_input(4);
		}

		// function for evaluation at single ip?
		inline void evaluate (number& value,
	        				const MathVector<dim>& globIP,
	        				number time, int si) const
		{
			UG_LOG("BinghamViscosityLinker::evaluate single called");
			number density;
			number viscosity;
			number yieldStress;
			MathMatrix<dim,dim> velocityGrad;

			(*m_spDensity)(density, globIP, time, si);
			(*m_spViscosity)(viscosity, globIP, time, si);
			(*m_spYieldStress)(yieldStress, globIP, time, si);
			(*m_spVelocityGrad)(velocityGrad, globIP, time, si);

			number innerSum = 0.0;

			// compute inner sum 
			for(int d1 = 0; d1 < dim; ++d1)
			{
				for(int d2 = 0; d2 < dim; ++d2)
				{
					innerSum += pow(velocityGrad(d1,d2) + velocityGrad(d2,d1),2);
				}
			}
			
			// compute mu = eta + sigma / \sqrt( delta + 1 / I)
			value = viscosity + yieldStress/sqrt(0.5+(1.0/(pow(2, dim))*innerSum));
			value *= 1./density;
		}

		// function for evaluation at multiple ips??
		template <int refDim>
		inline void evaluate(number vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GridObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			UG_LOG("BinghamViscosityLinker::evaluate called");
			std::vector<number> vDensity(nip);
			std::vector<number> vViscosity(nip);
			std::vector<number> vYieldStress(nip);
			std::vector<MathMatrix<dim,dim> > vVelocityGrad(nip);

			(*m_spDensity)(&vDensity[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spViscosity)(&vViscosity[0], vGlobIP, time, si,
								elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spYieldStress)(&vYieldStress[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spVelocityGrad)(&vVelocityGrad[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);

			for(size_t ip = 0; ip < nip; ++ip)
			{
				number mu = 0.0;

				// compute inner sum 
				for(int d1 = 0; d1 < dim; ++d1)
				{
					for(int d2 = 0; d2 < dim; ++d2)
					{
						mu += pow(vVelocityGrad[ip](d1,d2) + vVelocityGrad[ip](d2,d1),2);
					}
				}
				
				// compute mu = eta + sigma
				mu = vViscosity[ip] + vYieldStress[ip]/sqrt(0.5+(1.0/(pow(2, dim))*mu));
				vValue[ip] = mu * 1./vDensity[ip];
			}
		}

		template <int refDim>
		void eval_and_deriv(number vValue[],
				             const MathVector<dim> vGlobIP[],
				             number time, int si,
				             GridObject* elem,
				             const MathVector<dim> vCornerCoords[],
				             const MathVector<refDim> vLocIP[],
				             const size_t nip,
				             LocalVector* u,
				             bool bDeriv,
				             int s,
				             std::vector<std::vector<number> > vvvDeriv[],
				             const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		
		UG_LOG("BinghamViscosityLinker::eval_and_deriv called");
		//	get the data of the ip series
			const number* vDensity = m_spDensity->values(s);
			const number* vViscosity = m_spViscosity->values(s);
			const number* vYieldStress = m_spYieldStress->values(s);
			const MathMatrix<dim,dim>* vVelocityGrad = m_spVelocityGrad->values(s);			

			for(size_t ip = 0; ip < nip; ++ip)
			{
				number mu = 0.0;

				// compute mu = 
				for(int d1 = 0; d1 < dim; ++d1)
				{
					for(int d2 = 0; d2 < dim; ++d2)
					{
						mu += pow(vVelocityGrad[ip](d1,d2) + vVelocityGrad[ip](d2,d1),2);
					}
				}
				
				// compute mu = (eta + tau_F / \sqrt(delta + TrD^2)) / rho
				mu = vViscosity[ip] + vYieldStress[ip]/sqrt(0.5+(1.0/(pow(2, dim))*mu));
				vValue[ip] = mu * 1./vDensity[ip];
			}

			//	Compute the derivatives at all ips     //
			/////////////////////////////////////////////

		//	check if something is left to do
			if(!bDeriv || this->zero_derivative()) return;

		//	clear all derivative values
			this->set_zero(vvvDeriv, nip);

		//  Derivatives of Density
			if(m_spDDensity.valid() && !m_spDDensity->zero_derivative())
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
					for(size_t fct = 0; fct < m_spDDensity->num_fct(); ++fct)
					{
					//	get derivative of density w.r.t. to all functions
						const number* vDDensity = m_spDDensity->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = this->input_common_fct(_RHO_, fct);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
						{
							vvvDeriv[ip][commonFct][sh] -= vDDensity[sh] / vDensity[ip] * vValue[ip];
						}
					}
				}
			}

		//	Derivatives of Viscosity
			if(m_spDViscosity.valid() && !m_spDViscosity->zero_derivative())
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
					for(size_t fct = 0; fct < m_spDViscosity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const number* vDViscosity = m_spDViscosity->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = this->input_common_fct(_ETA_, fct);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
						{
							vvvDeriv[ip][commonFct][sh] += vDViscosity[sh] / vDensity[ip];
						}
					}
				}
			}

		//  Derivatives of yield stress
			if(m_spDYieldStress.valid() && !m_spDYieldStress->zero_derivative())
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
				//  Precompute 1 / (rho * \sqrt(delta + TrD^2))
					number rInvariant = vValue[ip] - (vViscosity[ip] / vDensity[ip]);

					for(size_t fct = 0; fct < m_spDYieldStress->num_fct(); ++fct)
					{
					//	get derivative of yield stress w.r.t. to all functions
						const number* vDYieldStress = m_spDYieldStress->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = this->input_common_fct(_TAU_, fct);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
						{
							vvvDeriv[ip][commonFct][sh] += vDYieldStress[sh] * rInvariant;
						}
					}
				}
			}

		//  Derivatives of velocity gradient
			UG_LOG("Derivatives of velocity gradient missing");
			/*if(m_spDVelocityGrad.valid() && !m_spDVelocityGrad->zero_derivative())
			{
				for(size_t ip = 0; ip < nip; ++ip)
				{
					for(size_t fct = 0; fct < m_spDYieldStress->num_fct(); ++fct)
					{
					//	get derivative of velocity gradient w.r.t. to all functions
						const number* vDYieldStress = m_spDYieldStress->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = this->input_common_fct(_DV_, fct);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < this->num_sh(commonFct); ++sh)
						{
							VecScaleAppend(vvvDeriv[ip][commonFct][sh], vDYieldStress[sh], rInvariant);
						}
					}
				}
			}*/
		}



	// Setter functions for imports
	///	set density import
		void set_density(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spDensity = data;
			m_spDDensity = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_RHO_, data, data);
		}

		void set_density(number val)
		{
			set_density(make_sp(new ConstUserNumber<dim>(val)));
		}

	///	set viscosity import
		void set_viscosity(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spViscosity = data;
			m_spDViscosity = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_ETA_, data, data);
		}

		void set_viscosity(number val)
		{
			set_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}

	///	set density import
		void set_yield_stress(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spYieldStress = data;
			m_spDYieldStress = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_TAU_, data, data);
		}

		void set_yield_stress(number val)
		{
			set_yield_stress(make_sp(new ConstUserNumber<dim>(val)));
		}

	/// set velocity gradient import
		void set_velocity_gradient(SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > data)
		{
			UG_LOG("Debug mark 1\n");
			m_spVelocityGrad = data;
			UG_LOG("Debug mark 2\n");
			try{
				UG_LOG(typeid(data).name() << "\n");
				m_spDVelocityGrad = data.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> >();
			} catch (std::exception& e) {
				UG_LOG("Error: " << e.what());
			}
			UG_LOG("Debug mark 3\n");
			base_type::set_input(_DV_, data, data);
		}

	protected:
		//  variables for storing imports
		///	import for density
			static const size_t _RHO_ = 0;
			SmartPtr<CplUserData<number, dim> > m_spDensity;
			SmartPtr<DependentUserData<number, dim> > m_spDDensity;
		///	import for viscosity
			static const size_t _ETA_ = 1;
			SmartPtr<CplUserData<number, dim> > m_spViscosity;
			SmartPtr<DependentUserData<number, dim> > m_spDViscosity;
		///	import for yield stress
			static const size_t _TAU_ = 2;
			SmartPtr<CplUserData<number, dim> > m_spYieldStress;
			SmartPtr<DependentUserData<number, dim> > m_spDYieldStress;
		///	import for velocity gradient
			static const size_t _DV_ = 3;
			SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > m_spVelocityGrad;
			SmartPtr<DependentUserData<MathMatrix<dim,dim>, dim> > m_spDVelocityGrad;
};

} // end of namespace ug

#endif // __H__UG__LIB_DISC__SPATIAL_DISC__BINGHAM_VISCOSITY_LINKER__

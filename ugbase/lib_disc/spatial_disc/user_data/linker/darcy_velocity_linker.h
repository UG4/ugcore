/*
 * darcy_velocity_linker.h
 *
 *      Author: andreasvogel
 */

#ifndef __H__UG__LIB_DISC__SPATIAL_DISC__DARCY_VELOCITY_LINKER__
#define __H__UG__LIB_DISC__SPATIAL_DISC__DARCY_VELOCITY_LINKER__

#include "linker.h"

namespace ug{


////////////////////////////////////////////////////////////////////////////////
// Darcy Velocity linker
////////////////////////////////////////////////////////////////////////////////

/// Hard Coded Linker for d3f
template <int dim>
class DarcyVelocityLinker
	: public StdDataLinker< DarcyVelocityLinker<dim>, MathVector<dim>, dim>
{
	///	Base class type
		typedef StdDataLinker< DarcyVelocityLinker<dim>, MathVector<dim>, dim> base_type;

	public:
		DarcyVelocityLinker() :
			m_spPermeability(NULL), m_spDPermeability(NULL),
			m_spViscosity(NULL), m_spDViscosity(NULL),
			m_spDensity(NULL), m_spDDensity(NULL),
			m_spGravity(NULL), m_spDGravity(NULL),
			m_spPressureGrad(NULL), m_spDPressureGrad(NULL)
		{
		//	this linker needs exactly five input
			this->set_num_input(5);
		}


		inline void evaluate (MathVector<dim>& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si) const
		{
			number density;
			number viscosity;
			MathVector<dim> gravity;
			MathVector<dim> pressureGrad;
			MathMatrix<dim,dim> permeability;

			(*m_spDensity)(density, globIP, time, si);
			(*m_spViscosity)(viscosity, globIP, time, si);
			(*m_spGravity)(gravity, globIP, time, si);
			(*m_spPressureGrad)(pressureGrad, globIP, time, si);
			(*m_spPermeability)(permeability, globIP, time, si);

		//	Variables
			MathVector<dim> Vel;

		//	compute rho*g
			VecScale(Vel, gravity, density);

		// 	compute rho*g - \nabla p
			VecSubtract(Vel, Vel, pressureGrad);

		//	compute Darcy velocity q := K / mu * (rho*g - \nabla p)
			MatVecMult(value, permeability, Vel);
			VecScale(value, value, 1./viscosity);
		}

		template <int refDim>
		inline void evaluate(MathVector<dim> vValue[],
		                     const MathVector<dim> vGlobIP[],
		                     number time, int si,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     LocalVector* u,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			std::vector<number> vDensity(nip);
			std::vector<number> vViscosity(nip);
			std::vector<MathVector<dim> > vGravity(nip);
			std::vector<MathVector<dim> > vPressureGrad(nip);
			std::vector<MathMatrix<dim,dim> > vPermeability(nip);

			(*m_spDensity)(&vDensity[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spViscosity)(&vViscosity[0], vGlobIP, time, si,
								elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spGravity)(&vGravity[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spPressureGrad)(&vPressureGrad[0], vGlobIP, time, si,
								elem, vCornerCoords, vLocIP, nip, u, vJT);
			(*m_spPermeability)(&vPermeability[0], vGlobIP, time, si,
							elem, vCornerCoords, vLocIP, nip, u, vJT);

			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	Variables
				MathVector<dim> Vel;

			//	compute rho*g
				VecScale(Vel, vGravity[ip], vDensity[ip]);

			// 	compute rho*g - \nabla p
				VecSubtract(Vel, Vel, vPressureGrad[ip]);

			//	compute Darcy velocity q := K / mu * (rho*g - \nabla p)
				MatVecMult(vValue[ip], vPermeability[ip], Vel);
				VecScale(vValue[ip], vValue[ip], 1./vViscosity[ip]);
			}
		}

		template <int refDim>
		void eval_and_deriv(MathVector<dim> vDarcyVel[],
		                    const MathVector<dim> vGlobIP[],
		                    number time, int si,
		                    GeometricObject* elem,
		                    const MathVector<dim> vCornerCoords[],
		                    const MathVector<refDim> vLocIP[],
		                    const size_t nip,
		                    LocalVector* u,
		                    bool bDeriv,
		                    int s,
		                    std::vector<std::vector<MathVector<dim> > > vvvDeriv[],
		                    const MathMatrix<refDim, dim>* vJT = NULL) const
		{
		//	get the data of the ip series
			const number* vDensity = m_spDensity->values(s);
			const number* vViscosity = m_spViscosity->values(s);
			const MathVector<dim>* vGravity = m_spGravity->values(s);
			const MathVector<dim>* vPressureGrad = m_spPressureGrad->values(s);
			const MathMatrix<dim,dim>* vPermeability = m_spPermeability->values(s);

			for(size_t ip = 0; ip < nip; ++ip)
			{
			//	Variables
				MathVector<dim> Vel;

			//	compute rho*g
				VecScale(Vel, vGravity[ip], vDensity[ip]);

			// 	compute rho*g - \nabla p
				VecSubtract(Vel, Vel, vPressureGrad[ip]);

			//	compute Darcy velocity q := K / mu * (rho*g - \nabla p)
				MatVecMult(vDarcyVel[ip], vPermeability[ip], Vel);
				VecScale(vDarcyVel[ip], vDarcyVel[ip], 1./vViscosity[ip]);
			}

			//	Compute the derivatives at all ips     //
			/////////////////////////////////////////////

		//	check if something to do
			if(!bDeriv || this->zero_derivative()) return;

		//	clear all derivative values
			this->set_zero(vvvDeriv, nip);

		//	Derivatives of Viscosity
			if(m_spDViscosity.valid() && !m_spDViscosity->zero_derivative())
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < m_spDViscosity->num_fct(); ++fct)
				{
				//	get derivative of viscosity w.r.t. to all functions
					const number* vDViscosity = m_spDViscosity->deriv(s, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_MU_, fct);

				//	loop all shapes and set the derivative
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
					//  DarcyVel_fct[sh] -= mu_fct_sh / mu * q
						VecScaleAppend(vvvDeriv[ip][commonFct][sh], -vDViscosity[sh] / vViscosity[ip], vDarcyVel[ip]);
					}
				}

		//	Derivatives of Density
			if(m_spDDensity.valid() && !m_spDDensity->zero_derivative())
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < m_spDDensity->num_fct(); ++fct)
				{
				//	get derivative of viscosity w.r.t. to all functions
					const number* vDDensity = m_spDDensity->deriv(s, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_RHO_, fct);

				//	Precompute K/mu * g
					MathVector<dim> Kmug;

				//	a) compute K * g
					MatVecMult(Kmug, vPermeability[ip], vGravity[ip]);

				//	b) compute K* g / mu
					VecScale(Kmug, Kmug, 1./vViscosity[ip]);

				//	loop all shapes and set the derivative
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
					//  DarcyVel_fct[sh] += K/mu * (rho_fct_sh * g)
						VecScaleAppend(vvvDeriv[ip][commonFct][sh],
									   vDDensity[sh], Kmug);
					}
				}

		//	Derivatives of Gravity
			if(m_spDGravity.valid() && !m_spDGravity->zero_derivative())
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < m_spDGravity->num_fct(); ++fct)
				{
				//	get derivative of viscosity w.r.t. to all functions
					const MathVector<dim>* vDGravity = m_spDGravity->deriv(s, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_G_, fct);

				//	Precompute K/mu * rho
					MathMatrix<dim,dim> Kmurho;

				//	a) compute K/mu * rho
					MatScale(Kmurho, vDensity[ip]/vViscosity[ip],vPermeability[ip]);

				//	loop all shapes and set the derivative
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						MathVector<dim> tmp;
						MatVecMult(tmp, Kmurho, vDGravity[sh]);

						vvvDeriv[ip][commonFct][sh] += tmp;
					}
				}

		//	Derivatives of Pressure
			if(m_spDPressureGrad.valid() && !m_spDPressureGrad->zero_derivative())
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < m_spDPressureGrad->num_fct(); ++fct)
				{
				//	get derivative of viscosity w.r.t. to all functions
					const MathVector<dim>* vDPressureGrad = m_spDPressureGrad->deriv(s, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_DP_, fct);

				//	Precompute -K/mu
					MathMatrix<dim,dim> Kmu;

				//	a) compute -K/mu
					MatScale(Kmu, -1.0/vViscosity[ip],vPermeability[ip]);

				//	loop all shapes and set the derivative
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						MathVector<dim> tmp;
						MatVecMult(tmp, Kmu, vDPressureGrad[sh]);

						vvvDeriv[ip][commonFct][sh] += tmp;
					}
				}

		//	Derivatives of Permeability
			if(m_spDPermeability.valid() && !m_spDPermeability->zero_derivative())
			for(size_t ip = 0; ip < nip; ++ip)
				for(size_t fct = 0; fct < m_spDPermeability->num_fct(); ++fct)
				{
				//	get derivative of viscosity w.r.t. to all functions
					const MathMatrix<dim,dim>* vDPermeability = m_spDPermeability->deriv(s, ip, fct);

				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(_K_, fct);

				//	Variables
					MathVector<dim> Vel;

				//	compute rho*g
					VecScale(Vel, vGravity[ip], vDensity[ip]);

				// 	compute rho*g - \nabla p
					VecSubtract(Vel, Vel, vPressureGrad[ip]);

				//	compute Darcy velocity q := K / mu * (rho*g - \nabla p)
					VecScale(Vel, Vel, 1./vViscosity[ip]);

				//	loop all shapes and set the derivative
					for(size_t sh = 0; sh < this->num_sh(fct); ++sh)
					{
						MathVector<dim> tmp;
						MatVecMult(tmp, vDPermeability[sh], Vel);

						vvvDeriv[ip][commonFct][sh] += tmp;
					}
				}
		}

	public:
	///	set permeability import
		void set_permeability(SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > data)
		{
			m_spPermeability = data;
			m_spDPermeability = data.template cast_dynamic<DependentUserData<MathMatrix<dim,dim>, dim> >();
			base_type::set_input(_K_, data, data);
		}

		void set_permeability(number val)
		{
			set_permeability(make_sp(new ConstUserMatrix<dim>(val)));
		}

	///	set permeability import
		void set_viscosity(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spViscosity = data;
			m_spDViscosity = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_MU_, data, data);
		}

		void set_viscosity(number val)
		{
			set_viscosity(make_sp(new ConstUserNumber<dim>(val)));
		}

	///	set density import
		void set_density(SmartPtr<CplUserData<number, dim> > data)
		{
			m_spDensity = data;
			m_spDDensity = data.template cast_dynamic<DependentUserData<number, dim> >();
			base_type::set_input(_RHO_, data, data);
		}

	///	set gravity import
		void set_gravity(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
		{
			m_spGravity = data;
			m_spDGravity = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
			base_type::set_input(_G_, data, data);
		}

	///	set pressure gradient import
		void set_pressure_gradient(SmartPtr<CplUserData<MathVector<dim>, dim> > data)
		{
			m_spPressureGrad = data;
			m_spDPressureGrad = data.template cast_dynamic<DependentUserData<MathVector<dim>, dim> >();
			base_type::set_input(_DP_, data, data);
		}

	protected:
	///	import for permeability
		static const size_t _K_ = 0;
		SmartPtr<CplUserData<MathMatrix<dim,dim>, dim> > m_spPermeability;
		SmartPtr<DependentUserData<MathMatrix<dim,dim>, dim> > m_spDPermeability;

	///	import for viscosity
		static const size_t _MU_ = 1;
		SmartPtr<CplUserData<number, dim> > m_spViscosity;
		SmartPtr<DependentUserData<number, dim> > m_spDViscosity;

	///	import for density
		static const size_t _RHO_ = 2;
		SmartPtr<CplUserData<number, dim> > m_spDensity;
		SmartPtr<DependentUserData<number, dim> > m_spDDensity;

	///	import for gravity
		static const size_t _G_ = 3;
		SmartPtr<CplUserData<MathVector<dim>, dim> > m_spGravity;
		SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDGravity;

	///	import for pressure gradient
		static const size_t _DP_ = 4;
		SmartPtr<CplUserData<MathVector<dim>, dim> > m_spPressureGrad;
		SmartPtr<DependentUserData<MathVector<dim>, dim> > m_spDPressureGrad;
};

} // end namespace ug

#endif /* __H__UG__LIB_DISC__SPATIAL_DISC__DARCY_VELOCITY_LINKER__ */

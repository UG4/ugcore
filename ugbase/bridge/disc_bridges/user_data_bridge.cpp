// author: Andreas Vogel

#include <iostream>
#include <sstream>
#include <boost/function.hpp>

// include bridge
#include "bridge/bridge.h"
#include "bridge/util.h"

#include "common/common.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_linker.h"
#include "lib_disc/spatial_disc/user_data/user_function.h"
#include "lib_disc/spatial_disc/user_data/scale_add_linker.h"
#include "lib_disc/spatial_disc/user_data/inverse_linker.h"

using namespace std;

namespace ug{
namespace bridge{

/// Hard Coded Linker for d3f
template <int dim>
class DarcyVelocityLinker
	: public StdDataLinker< DarcyVelocityLinker<dim>, MathVector<dim>, dim>
{
	///	Base class type
		typedef DataLinker<MathVector<dim>, dim> base_type;

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
		inline void evaluate (MathVector<dim>& value,
		                      const MathVector<dim>& globIP,
		                      number time, int si,
		                      LocalVector& u,
		                      GeometricObject* elem,
		                      const MathVector<dim> vCornerCoords[],
		                      const MathVector<refDim>& locIP) const
		{
			number density;
			number viscosity;
			MathVector<dim> gravity;
			MathVector<dim> pressureGrad;
			MathMatrix<dim,dim> permeability;

			(*m_spDensity)(density, globIP, time, si, u,
							elem, vCornerCoords, locIP);
			(*m_spViscosity)(viscosity, globIP, time, si, u,
								elem, vCornerCoords, locIP);
			(*m_spGravity)(gravity, globIP, time, si, u,
							elem, vCornerCoords, locIP);
			(*m_spPressureGrad)(pressureGrad, globIP, time, si, u,
								elem, vCornerCoords, locIP);
			(*m_spPermeability)(permeability, globIP, time, si, u,
							elem, vCornerCoords, locIP);

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
		                     LocalVector& u,
		                     GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[],
		                     const MathVector<refDim> vLocIP[],
		                     const size_t nip,
		                     const MathMatrix<refDim, dim>* vJT = NULL) const
		{
			std::vector<number> vDensity(nip);
			std::vector<number> vViscosity(nip);
			std::vector<MathVector<dim> > vGravity(nip);
			std::vector<MathVector<dim> > vPressureGrad(nip);
			std::vector<MathMatrix<dim,dim> > vPermeability(nip);

			(*m_spDensity)(&vDensity[0], vGlobIP, time, si, u,
							elem, vCornerCoords, vLocIP, nip, vJT);
			(*m_spViscosity)(&vViscosity[0], vGlobIP, time, si, u,
								elem, vCornerCoords, vLocIP, nip, vJT);
			(*m_spGravity)(&vGravity[0], vGlobIP, time, si, u,
							elem, vCornerCoords, vLocIP, nip, vJT);
			(*m_spPressureGrad)(&vPressureGrad[0], vGlobIP, time, si, u,
								elem, vCornerCoords, vLocIP, nip, vJT);
			(*m_spPermeability)(&vPermeability[0], vGlobIP, time, si, u,
							elem, vCornerCoords, vLocIP, nip, vJT);

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

		virtual void compute(LocalVector* u, GeometricObject* elem,
		                     const MathVector<dim> vCornerCoords[], bool bDeriv = false)
		{
		//	Compute the Darcy velocity at all ips  //
		/////////////////////////////////////////////

		//	loop all time series and every integration point of the series
			for(size_t s = 0; s < this->num_series(); ++s)
			{
			//	get the data of the ip series
				const number* vDensity = m_spDensity->values(s);
				const number* vViscosity = m_spViscosity->values(s);
				const MathVector<dim>* vGravity = m_spGravity->values(s);
				const MathVector<dim>* vPressureGrad = m_spPressureGrad->values(s);
				const MathMatrix<dim,dim>* vPermeability = m_spPermeability->values(s);

			//	get the data to be filled
				MathVector<dim>* vDarcyVel = this->values(s);

				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
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
			}

			//	Compute the derivatives at all ips     //
			/////////////////////////////////////////////

		//	check if something to do
			if(!bDeriv || this->zero_derivative()) return;

		//	clear all derivative values
			this->clear_derivative_values();

		//	loop all series
			for(size_t s = 0; s < this->num_series(); ++s)
			{
			//	get the data of the ip series
				const number* vDensity = m_spDensity->values(s);
				const number* vViscosity = m_spViscosity->values(s);
				const MathVector<dim>* vGravity = m_spGravity->values(s);
				const MathVector<dim>* vPressureGrad = m_spPressureGrad->values(s);
				const MathMatrix<dim,dim>* vPermeability = m_spPermeability->values(s);

			//	get the data to be filled
				MathVector<dim>* vDarcyVel = this->values(s);

			//	Derivatives of Viscosity
				if(m_spDViscosity.valid() && !m_spDViscosity->zero_derivative())
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
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
							VecSubtract(this->deriv(s, ip, commonFct, sh), vDarcyVel[ip], -vDViscosity[sh] / vViscosity[ip]);
						}
					}

			//	Derivatives of Density
				if(m_spDDensity.valid() && !m_spDDensity->zero_derivative())
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
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
							VecScaleAppend(this->deriv(s, ip, commonFct, sh),
							               vDDensity[sh], Kmug);
						}
					}

			//	Derivatives of Gravity
				if(m_spDGravity.valid() && !m_spDGravity->zero_derivative())
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
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

							this->deriv(s, ip, commonFct, sh) += tmp;
						}
					}

			//	Derivatives of Pressure
				if(m_spDPressureGrad.valid() && !m_spDPressureGrad->zero_derivative())
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
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

							this->deriv(s, ip, commonFct, sh) += tmp;
						}
					}

			//	Derivatives of Permeability
				if(m_spDPermeability.valid() && !m_spDPermeability->zero_derivative())
				for(size_t ip = 0; ip < this->num_ip(s); ++ip)
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

							this->deriv(s, ip, commonFct, sh) += tmp;
						}
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
			set_permeability(CreateSmartPtr(new ConstUserMatrix<dim>(val)));
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
			set_viscosity(CreateSmartPtr(new ConstUserNumber<dim>(val)));
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


template <typename TData, int dim>
void RegisterUserDataType(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	string type = user_data_traits<TData>::name();

//	User"Type"
//	NOTE: For better readability this class is named User"Type"
//	      in vrl and lua. E.g. UserNumber, UserVector, ...
	{
		typedef UserData<TData, dim> T;
		typedef UserDataInfo TBase1;
		string name = string("User").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp)
			.add_method("get_dim", &T::get_dim)
			.add_method("type", &T::type);
		reg.add_class_to_group(name, string("User").append(type), dimTag);
	}

//	CondUser"Type"
//	NOTE: For better readability this class is named CondUser"Type"
//	 	  in vrl and lua. E.g. CondUserNumber, CondUserVector, ...
	{
		typedef UserData<TData, dim, bool> T;
		typedef UserDataInfo TBase1;
		string name = string("CondUser").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp);
		reg.add_class_to_group(name, string("CondUser").append(type), dimTag);
	}

//	CplUser"Type"
	{
		typedef CplUserData<TData, dim> T;
		typedef UserData<TData,dim> TBase1;
		string name = string("CplUser").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp)
			.add_method("get_dim", &T::get_dim)
			.add_method("type", &T::type);
		reg.add_class_to_group(name, string("CplUser").append(type), dimTag);
	}

//	CondCplUser"Type"
	{
		typedef CplUserData<TData, dim, bool> T;
		typedef UserData<TData,dim,bool> TBase1;
		string name = string("CondCplUser").append(type).append(dimSuffix);
		reg.add_class_<T,TBase1>(name, grp);
		reg.add_class_to_group(name, string("CondCplUser").append(type), dimTag);
	}

//	DependentUserData"Type"
	{
		typedef DependentUserData<TData, dim> T;
		typedef CplUserData<TData, dim> TBase;
		string name = string("DependentUserData").append(type).append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, string("DependentUserData").append(type), dimTag);
	}

//	DataLinker"Type"
	{
		typedef DataLinker<TData, dim> T;
		typedef DependentUserData<TData, dim> TBase;
		string name = string("DataLinker").append(type).append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, string("DataLinker").append(type), dimTag);
	}

//	ScaleAddLinker"Type"
	{
		typedef ScaleAddLinker<TData, dim, number> T;
		typedef DataLinker<TData, dim> TBase;
		string name = string("ScaleAddLinker").append(type).append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number , SmartPtr<CplUserData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , number)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number,number)>(&T::add))
			.add_constructor()
			.template add_constructor<void (*)(const ScaleAddLinker<TData, dim, number>&)>()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("ScaleAddLinker").append(type), dimTag);
	}



}

namespace UserDataBridge{

/**
 * Class exporting the functionality. All functionality that is to
 * be used in scripts or visualization must be registered here.
 */
struct Functionality
{

/**
 * Function called for the registration of Dimension dependent parts.
 * All Functions and Classes depending on the Dimension
 * are to be placed here when registering. The method is called for all
 * available Dimension types, based on the current build options.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
template <int dim>
static void Dimension(Registry& reg, string grp)
{
	string dimSuffix = GetDimensionSuffix<dim>();
	string dimTag = GetDimensionTag<dim>();

	RegisterUserDataType<number, dim>(reg, grp);
	RegisterUserDataType<MathVector<dim>, dim>(reg, grp);
	RegisterUserDataType<MathMatrix<dim,dim>, dim>(reg, grp);
	RegisterUserDataType<MathTensor<4,dim>, dim>(reg, grp);

//	ConstUserNumber
	{
		typedef ConstUserNumber<dim> T;
		typedef CplUserData<number, dim> TBase;
		string name = string("ConstUserNumber").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(number)>("Value")
			.add_method("set", &T::set, "", "Value")
			.add_method("print", &T::print)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstUserNumber", dimTag);
	}

//	ConstUserVector
	{
		typedef ConstUserVector<dim> T;
		typedef CplUserData<MathVector<dim>, dim> TBase;
		string name = string("ConstUserVector").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(number)>("Values")
			.add_method("set_all_entries", &T::set_all_entries)
			.add_method("set_entry", &T::set_entry)
			.add_method("print", &T::print)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstUserVector", dimTag);
	}

//	ConstUserMatrix
	{
		typedef ConstUserMatrix<dim> T;
		typedef CplUserData<MathMatrix<dim, dim>, dim> TBase;
		string name = string("ConstUserMatrix").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.template add_constructor<void (*)(number)>("Diagonal Value")
			.add_method("set_diag_tensor", &T::set_diag_tensor)
			.add_method("set_all_entries", &T::set_all_entries)
			.add_method("set_entry", &T::set_entry)
			.add_method("print", &T::print)
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ConstUserMatrix", dimTag);
	}

//	DarcyVelocityLinker
	{
		typedef DarcyVelocityLinker<dim> T;
		typedef DataLinker<MathVector<dim>, dim> TBase;
		string name = string("DarcyVelocityLinker").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("set_density", &T::set_density)
			.add_method("set_gravity", &T::set_gravity)
			.add_method("set_permeability", static_cast<void (T::*)(number)>(&T::set_permeability))
			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<CplUserData<MathMatrix<dim,dim>,dim> >)>(&T::set_permeability))
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity))
			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> >)>(&T::set_viscosity))
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DarcyVelocityLinker", dimTag);
	}
	//	InverseLinker"Type"
		{
			typedef InverseLinker<dim> T;
			typedef DataLinker<number,dim> TBase;

				string name = string("InverseLinker").append(dimSuffix);
				reg.add_class_<T,TBase>(name, grp)
				.add_method("divide", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , SmartPtr<CplUserData<number,dim> >)>(&T::divide))
				.add_method("divide", static_cast<void (T::*)(number , SmartPtr<CplUserData<number,dim> >)>(&T::divide))
				.add_method("divide", static_cast<void (T::*)(SmartPtr<CplUserData<number,dim> > , number)>(&T::divide))
				.add_method("divide", static_cast<void (T::*)(number,number)>(&T::divide))
				.add_constructor()
				.template add_constructor<void (*)(const InverseLinker<dim>&)>()
				.set_construct_as_smart_pointer(true);
			reg.add_class_to_group(name, string("InverseLinker"), dimTag);

		}

}

/**
 * Function called for the registration of Domain and Algebra independent parts.
 * All Functions and Classes not depending on Domain and Algebra
 * are to be placed here when registering.
 *
 * @param reg				registry
 * @param parentGroup		group for sorting of functionality
 */
static void Common(Registry& reg, string grp)
{
	reg.add_class_<IFunction<number> >("IFunctionNumber", grp);

//	UserDataInfo
	{
		reg.add_class_<UserDataInfo>("UserDataInfo", grp);
	}

}

}; // end Functionality
}// end UserData

void RegisterBridge_UserData(Registry& reg, string grp)
{
//	get group string
	grp.append("/Discretization/SpatialDisc/UserData");
	typedef UserDataBridge::Functionality Functionality;

	try{
		RegisterCommon<Functionality>(reg,grp);
		RegisterDimensionDependent<Functionality>(reg,grp);
	}
	UG_REGISTRY_CATCH_THROW(grp);
}

} // end namespace
} // end namepace

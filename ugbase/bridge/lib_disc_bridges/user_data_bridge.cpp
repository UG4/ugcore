// author: Andreas Vogel

#include <iostream>
#include <sstream>
#include <boost/function.hpp>

#include "bridge/bridge.h"

#include "common/common.h"
#include "lib_disc/spatial_disc/ip_data/const_user_data.h"
#include "lib_disc/spatial_disc/ip_data/data_linker.h"
#include "lib_disc/spatial_disc/ip_data/user_function.h"

using namespace std;

namespace ug
{
namespace bridge
{

/// Hard Coded Linker for d3f
template <int dim>
class ElderDensityLinker
	: public DataLinkerEqualData<number, dim, number>
{
	///	Base class type
		typedef DataLinkerEqualData<number, dim, number> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::time;

	//	explicitly forward methods of IPData
		using base_type::ip;
		using base_type::value;

	//	explicitly forward methods of IDependentIPData
		using base_type::num_fct;

	//	explicitly forward methods of DependentIPData
		using base_type::num_sh;
		using base_type::deriv;

	//	explicitly forward methods of Data Linker
		using base_type::set_num_input;
		using base_type::input_value;
		using base_type::input_deriv;

	public:
		ElderDensityLinker()
		{
		//	this linker needs exactly one input
			set_num_input(1);
		}

		virtual void compute(bool bDeriv)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
					value(s, ip) = 1e3 + 0.2e3 * input_value(0, s, ip);
				}

			if(!bDeriv || this->zero_derivative()) return;

			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < num_fct(); ++fct)
						for(size_t dof = 0; dof < num_sh(fct); ++dof)
						{
							deriv(s, ip, fct, dof) = 0.2e3 * input_deriv(0, s, ip, fct, dof);
						}
		}
};


/// Hard Coded Linker for d3f
template <int dim>
class DarcyVelocityLinker
	: public DataLinker<MathVector<dim>, dim>
{
	///	Base class type
		typedef DataLinker<MathVector<dim>, dim> base_type;

	//	explicitly forward methods of IIPData
		using base_type::num_series;
		using base_type::num_ip;
		using base_type::time;

	//	explicitly forward methods of IPData
		using base_type::ip;
		using base_type::value;
		using base_type::values;

	//	explicitly forward methods of IDependentIPData
		using base_type::num_fct;

	//	explicitly forward methods of DependentIPData
		using base_type::num_sh;
		using base_type::deriv;

	//	explicitly forward methods of Data Linker
		using base_type::set_num_input;
		using base_type::input_common_fct;

	public:
		DarcyVelocityLinker() :
			m_spPermeability(NULL), m_spDPermeability(NULL),
			m_spViscosity(NULL), m_spDViscosity(NULL),
			m_spDensity(NULL), m_spDDensity(NULL),
			m_spGravity(NULL), m_spDGravity(NULL),
			m_spPressureGrad(NULL), m_spDPressureGrad(NULL)
		{
		//	this linker needs exactly one input
			set_num_input(5);
		}

		virtual void compute(bool bDeriv)
		{
		//	Compute the Darcy velocity at all ips  //
		/////////////////////////////////////////////

		//	loop all time series and every integration point of the series
			for(size_t s = 0; s < num_series(); ++s)
			{
			//	get the data of the ip series
				const number* vDensity = m_spDensity->values(s);
				const number* vViscosity = m_spViscosity->values(s);
				const MathVector<dim>* vGravity = m_spGravity->values(s);
				const MathVector<dim>* vPressureGrad = m_spPressureGrad->values(s);
				const MathMatrix<dim,dim>* vPermeability = m_spPermeability->values(s);

			//	get the data to be filled
				MathVector<dim>* vDarcyVel = values(s);

				for(size_t ip = 0; ip < num_ip(s); ++ip)
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
			for(size_t s = 0; s < num_series(); ++s)
			{
			//	get the data of the ip series
				const number* vDensity = m_spDensity->values(s);
				const number* vViscosity = m_spViscosity->values(s);
				const MathVector<dim>* vGravity = m_spGravity->values(s);
				const MathVector<dim>* vPressureGrad = m_spPressureGrad->values(s);
				const MathMatrix<dim,dim>* vPermeability = m_spPermeability->values(s);

			//	get the data to be filled
				MathVector<dim>* vDarcyVel = values(s);

			//	Derivatives of Viscosity
				if(m_spDViscosity.valid() && !m_spDViscosity->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_spDViscosity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const number* vDViscosity = m_spDViscosity->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = input_common_fct(_MU_, fct);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < num_sh(fct); ++sh)
						{
						//  DarcyVel_fct[sh] -= mu_fct_sh / mu * q
							VecSubtract(deriv(s, ip, commonFct, sh), vDarcyVel[ip], -vDViscosity[sh] / vViscosity[ip]);
						}
					}

			//	Derivatives of Density
				if(m_spDDensity.valid() && !m_spDDensity->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_spDDensity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const number* vDDensity = m_spDDensity->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = input_common_fct(_RHO_, fct);

					//	Precompute K/mu * g
						MathVector<dim> Kmug;

					//	a) compute K * g
						MatVecMult(Kmug, vPermeability[ip], vGravity[ip]);

					//	b) compute K* g / mu
						VecScale(Kmug, Kmug, 1./vViscosity[ip]);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < num_sh(fct); ++sh)
						{
						//  DarcyVel_fct[sh] += K/mu * (rho_fct_sh * g)
							VecScaleAppend(deriv(s, ip, commonFct, sh),
							               vDDensity[sh], Kmug);
						}
					}

			//	Derivatives of Gravity
				if(m_spDGravity.valid() && !m_spDGravity->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_spDGravity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const MathVector<dim>* vDGravity = m_spDGravity->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = input_common_fct(_G_, fct);

					//	Precompute K/mu * rho
						MathMatrix<dim,dim> Kmurho;

					//	a) compute K/mu * rho
						MatScale(Kmurho, vDensity[ip]/vViscosity[ip],vPermeability[ip]);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < num_sh(fct); ++sh)
						{
							MathVector<dim> tmp;
							MatVecMult(tmp, Kmurho, vDGravity[sh]);

							deriv(s, ip, commonFct, sh) += tmp;
						}
					}

			//	Derivatives of Pressure
				if(m_spDPressureGrad.valid() && !m_spDPressureGrad->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_spDPressureGrad->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const MathVector<dim>* vDPressureGrad = m_spDPressureGrad->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = input_common_fct(_DP_, fct);

					//	Precompute -K/mu
						MathMatrix<dim,dim> Kmu;

					//	a) compute -K/mu
						MatScale(Kmu, -1.0/vViscosity[ip],vPermeability[ip]);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < num_sh(fct); ++sh)
						{
							MathVector<dim> tmp;
							MatVecMult(tmp, Kmu, vDPressureGrad[sh]);

							deriv(s, ip, commonFct, sh) += tmp;
						}
					}

			//	Derivatives of Permeability
				if(m_spDPermeability.valid() && !m_spDPermeability->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_spDPermeability->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const MathMatrix<dim,dim>* vDPermeability = m_spDPermeability->deriv(s, ip, fct);

					//	get common fct id for this function
						const size_t commonFct = input_common_fct(_K_, fct);

					//	Variables
						MathVector<dim> Vel;

					//	compute rho*g
						VecScale(Vel, vGravity[ip], vDensity[ip]);

					// 	compute rho*g - \nabla p
						VecSubtract(Vel, Vel, vPressureGrad[ip]);

					//	compute Darcy velocity q := K / mu * (rho*g - \nabla p)
						VecScale(Vel, Vel, 1./vViscosity[ip]);

					//	loop all shapes and set the derivative
						for(size_t sh = 0; sh < num_sh(fct); ++sh)
						{
							MathVector<dim> tmp;
							MatVecMult(tmp, vDPermeability[sh], Vel);

							deriv(s, ip, commonFct, sh) += tmp;
						}
					}

			}
		}

	public:
	///	set permeability import
		void set_permeability(SmartPtr<IPData<MathMatrix<dim,dim>, dim> > data)
		{
			m_spPermeability = data;
			m_spDPermeability = data.template cast_dynamic<DependentIPData<MathMatrix<dim,dim>, dim> >();
			base_type::set_input(_K_, data);
		}

		void set_permeability(number val)
		{
			set_permeability(CreateSmartPtr(new ConstUserMatrix<dim>(val)));
		}

	///	set permeability import
		void set_viscosity(SmartPtr<IPData<number, dim> > data)
		{
			m_spViscosity = data;
			m_spDViscosity = data.template cast_dynamic<DependentIPData<number, dim> >();
			base_type::set_input(_MU_, data);
		}

		void set_viscosity(number val)
		{
			set_viscosity(CreateSmartPtr(new ConstUserNumber<dim>(val)));
		}

	///	set density import
		void set_density(SmartPtr<IPData<number, dim> > data)
		{
			m_spDensity = data;
			m_spDDensity = data.template cast_dynamic<DependentIPData<number, dim> >();
			base_type::set_input(_RHO_, data);
		}

	///	set gravity import
		void set_gravity(SmartPtr<IPData<MathVector<dim>, dim> > data)
		{
			m_spGravity = data;
			m_spDGravity = data.template cast_dynamic<DependentIPData<MathVector<dim>, dim> >();
			base_type::set_input(_G_, data);
		}

	///	set pressure gradient import
		void set_pressure_gradient(SmartPtr<IPData<MathVector<dim>, dim> > data)
		{
			m_spPressureGrad = data;
			m_spDPressureGrad = data.template cast_dynamic<DependentIPData<MathVector<dim>, dim> >();
			base_type::set_input(_DP_, data);
		}

	protected:
	///	import for permeability
		static const size_t _K_ = 0;
		SmartPtr<IPData<MathMatrix<dim,dim>, dim> > m_spPermeability;
		SmartPtr<DependentIPData<MathMatrix<dim,dim>, dim> > m_spDPermeability;

	///	import for viscosity
		static const size_t _MU_ = 1;
		SmartPtr<IPData<number, dim> > m_spViscosity;
		SmartPtr<DependentIPData<number, dim> > m_spDViscosity;

	///	import for density
		static const size_t _RHO_ = 2;
		SmartPtr<IPData<number, dim> > m_spDensity;
		SmartPtr<DependentIPData<number, dim> > m_spDDensity;

	///	import for gravity
		static const size_t _G_ = 3;
		SmartPtr<IPData<MathVector<dim>, dim> > m_spGravity;
		SmartPtr<DependentIPData<MathVector<dim>, dim> > m_spDGravity;

	///	import for pressure gradient
		static const size_t _DP_ = 4;
		SmartPtr<IPData<MathVector<dim>, dim> > m_spPressureGrad;
		SmartPtr<DependentIPData<MathVector<dim>, dim> > m_spDPressureGrad;
};


template <typename TData, int dim>
bool RegisterUserDataType(Registry& reg, string type, string parentGroup)
{
	string grp = string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	User"Type"
//	NOTE: For better readability this class is named User"Type"
//	      in vrl and lua. E.g. UserNumber, UserVector, ...
	{
		typedef IPData<TData, dim> T;
		string name = string("User").append(type).append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, string("User").append(type), dimTag);
	}

//	CondUser"Type"
//	NOTE: For better readability this class is named CondUser"Type"
//	 	  in vrl and lua. E.g. CondUserNumber, CondUserVector, ...
	{
		typedef IPData<TData, dim, bool> T;
		string name = string("CondUser").append(type).append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, string("CondUser").append(type), dimTag);
	}

//	DependentIPData"Type"
	{
		typedef DependentIPData<TData, dim> T;
		typedef IPData<TData, dim> TBase;
		string name = string("DependentIPData").append(type).append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, string("DependentIPData").append(type), dimTag);
	}

//	DataLinker"Type"
	{
		typedef DataLinker<TData, dim> T;
		typedef DependentIPData<TData, dim> TBase;
		string name = string("DataLinker").append(type).append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp);
		reg.add_class_to_group(name, string("DataLinker").append(type), dimTag);
	}

//	DataLinkerEqual"Type"
	{
		typedef DataLinkerEqualData<TData, dim, number> T;
		typedef DataLinker<TData, dim> TBase;
		string name = string("DataLinkerEqual").append(type).append(dimSuffix);
		reg.add_class_<T, TBase >(name, grp)
			.add_method("set_input", &T::set_input);
		reg.add_class_to_group(name, string("DataLinkerEqual").append(type), dimTag);
	}

//	ScaleAddLinker"Type"
	{
		typedef ScaleAddLinker<TData, dim, number> T;
		typedef DataLinker<TData, dim> TBase;
		string name = string("ScaleAddLinker").append(type).append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<number,dim> > , SmartPtr<IPData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number , SmartPtr<IPData<TData,dim> >)>(&T::add))
			.add_method("add", static_cast<void (T::*)(SmartPtr<IPData<number,dim> > , number)>(&T::add))
			.add_method("add", static_cast<void (T::*)(number,number)>(&T::add))
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, string("ScaleAddLinker").append(type), dimTag);
	}

	return true;
}

template <int dim>
bool RegisterUserData(Registry& reg, string parentGroup)
{
	string grp = std::string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

	RegisterUserDataType<number, dim>(reg, "Number", parentGroup);
	RegisterUserDataType<MathVector<dim>, dim>(reg, "Vector", parentGroup);
	RegisterUserDataType<MathMatrix<dim,dim>, dim>(reg, "Matrix", parentGroup);
	RegisterUserDataType<MathTensor<4,dim>, dim>(reg, "Tensor4", parentGroup);

//	ConstUserNumber
	{
		typedef ConstUserNumber<dim> T;
		typedef IPData<number, dim> TBase;
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
		typedef IPData<MathVector<dim>, dim> TBase;
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
		typedef IPData<MathMatrix<dim, dim>, dim> TBase;
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

//	ElderDensityLinker
	{
		typedef ElderDensityLinker<dim> T;
		typedef DataLinkerEqualData<number, dim, number> TBase;
		string name = string("ElderDensityLinker").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "ElderDensityLinker", dimTag);
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
			.add_method("set_permeability", static_cast<void (T::*)(SmartPtr<IPData<MathMatrix<dim,dim>,dim> >)>(&T::set_permeability))
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", static_cast<void (T::*)(number)>(&T::set_viscosity))
			.add_method("set_viscosity", static_cast<void (T::*)(SmartPtr<IPData<number,dim> >)>(&T::set_viscosity))
			.add_constructor()
			.set_construct_as_smart_pointer(true);
		reg.add_class_to_group(name, "DarcyVelocityLinker", dimTag);
	}

	return true;
}

bool RegisterUserData(Registry& reg, string parentGroup)
{
//	get group string
	std::string grp = parentGroup; grp.append("/Discretization/SpatialDisc/UserData");

#ifdef UG_DIM_1
	RegisterUserData<1>(reg, grp);
#endif
#ifdef UG_DIM_2
	RegisterUserData<2>(reg, grp);
#endif
#ifdef UG_DIM_3
	RegisterUserData<3>(reg, grp);
#endif

	reg.add_class_<IFunction<number> >("IFunctionNumber", grp);

	return true;
}

} // end namespace
} // end namepace

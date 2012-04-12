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
			m_pPermeability(NULL), m_pdPermeability(NULL),
			m_pViscosity(NULL), m_pdViscosity(NULL),
			m_pDensity(NULL), m_pdDensity(NULL),
			m_pGravity(NULL), m_pdGravity(NULL),
			m_pPressureGrad(NULL), m_pdPressureGrad(NULL)
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
				const number* vDensity = m_pDensity->values(s);
				const number* vViscosity = m_pViscosity->values(s);
				const MathVector<dim>* vGravity = m_pGravity->values(s);
				const MathVector<dim>* vPressureGrad = m_pPressureGrad->values(s);
				const MathMatrix<dim,dim>* vPermeability = m_pPermeability->values(s);

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
				const number* vDensity = m_pDensity->values(s);
				const number* vViscosity = m_pViscosity->values(s);
				const MathVector<dim>* vGravity = m_pGravity->values(s);
				const MathVector<dim>* vPressureGrad = m_pPressureGrad->values(s);
				const MathMatrix<dim,dim>* vPermeability = m_pPermeability->values(s);

			//	get the data to be filled
				MathVector<dim>* vDarcyVel = values(s);

			//	Derivatives of Viscosity
				if(m_pdViscosity != NULL && !m_pdViscosity->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_pdViscosity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const number* vDViscosity = m_pdViscosity->deriv(s, ip, fct);

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
				if(m_pdDensity != NULL && !m_pdDensity->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_pdDensity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const number* vDDensity = m_pdDensity->deriv(s, ip, fct);

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
				if(m_pdGravity != NULL && !m_pdGravity->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_pdGravity->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const MathVector<dim>* vDGravity = m_pdGravity->deriv(s, ip, fct);

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
				if(m_pdPressureGrad != NULL && !m_pdPressureGrad->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_pdPressureGrad->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const MathVector<dim>* vDPressureGrad = m_pdPressureGrad->deriv(s, ip, fct);

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
				if(m_pdPermeability != NULL && !m_pdPermeability->zero_derivative())
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < m_pdPermeability->num_fct(); ++fct)
					{
					//	get derivative of viscosity w.r.t. to all functions
						const MathMatrix<dim,dim>* vDPermeability = m_pdPermeability->deriv(s, ip, fct);

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
		void set_permeability(IPData<MathMatrix<dim,dim>, dim>& data)
		{
			m_pPermeability = &data;
			m_pdPermeability = dynamic_cast<DependentIPData<MathMatrix<dim,dim>, dim>*>(&data);
			base_type::set_input(_K_, &data);
		}

	///	set permeability import
		void set_viscosity(IPData<number, dim>& data)
		{
			m_pViscosity = &data;
			m_pdViscosity = dynamic_cast<DependentIPData<number, dim>*>(&data);
			base_type::set_input(_MU_, &data);
		}

	///	set density import
		void set_density(IPData<number, dim>& data)
		{
			m_pDensity = &data;
			m_pdDensity = dynamic_cast<DependentIPData<number, dim>*>(&data);
			base_type::set_input(_RHO_, &data);
		}

	///	set gravity import
		void set_gravity(IPData<MathVector<dim>, dim>& data)
		{
			m_pGravity = &data;
			m_pdGravity = dynamic_cast<DependentIPData<MathVector<dim>, dim>*>(&data);
			base_type::set_input(_G_, &data);
		}

	///	set pressure gradient import
		void set_pressure_gradient(IPData<MathVector<dim>, dim>& data)
		{
			m_pPressureGrad = &data;
			m_pdPressureGrad = dynamic_cast<DependentIPData<MathVector<dim>, dim>*>(&data);
			base_type::set_input(_DP_, &data);
		}

	protected:
	///	import for permeability
		static const size_t _K_ = 0;
		IPData<MathMatrix<dim,dim>, dim>* m_pPermeability;
		DependentIPData<MathMatrix<dim,dim>, dim>* m_pdPermeability;

	///	import for viscosity
		static const size_t _MU_ = 1;
		IPData<number, dim>* m_pViscosity;
		DependentIPData<number, dim>* m_pdViscosity;

	///	import for density
		static const size_t _RHO_ = 2;
		IPData<number, dim>* m_pDensity;
		DependentIPData<number, dim>* m_pdDensity;

	///	import for gravity
		static const size_t _G_ = 3;
		IPData<MathVector<dim>, dim>* m_pGravity;
		DependentIPData<MathVector<dim>, dim>* m_pdGravity;

	///	import for pressure gradient
		static const size_t _DP_ = 4;
		IPData<MathVector<dim>, dim>* m_pPressureGrad;
		DependentIPData<MathVector<dim>, dim>* m_pdPressureGrad;
};


template <typename TData, int dim>
bool RegisterUserDataType(Registry& reg, string type, string parentGroup)
{
	string grp = string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	IPData"Type"
//	User"Type" NOTE: For better readability this class is named User"Type"
//				     in vrl and lua. E.g. UserNumber, UserVector, ...
	{
		typedef IPData<TData, dim> T;
		string name = string("User").append(type).append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, string("User").append(type), dimTag);
	}

//	IPCondData"Type"
	{
		typedef IPData<TData, dim, bool> T;
		string name = string("IPCondData").append(type).append(dimSuffix);
		reg.add_class_<T>(name, grp);
		reg.add_class_to_group(name, string("IPCondData").append(type), dimTag);
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
			.add_method("add", &T::add)
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
			.add_method("set_permeability", &T::set_permeability)
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", &T::set_viscosity)
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

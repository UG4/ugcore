

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"
#include "lib_discretization/spatial_discretization/ip_data/data_linker.h"
#include "ug_script/ug_script.h"
#include <iostream>
#include <sstream>

namespace ug
{
namespace bridge
{

class PrintUserNumber2d
{
	protected:
		typedef IUserData<number, 2>::functor_type NumberFunctor;

	public:
		void set_user_number(IUserData<number, 2>& user)
		{
			m_Number = user.get_functor();
		}

		number print(number x, number y)
		{
			MathVector<2> v(x,y);
			number time = 0.0;
			number ret;

			if(m_Number)
				m_Number(ret, v, time);
			else
			{
				UG_LOG("Functor not set. \n");
				ret = -1;
			}

			return ret;
		}

	private:
		NumberFunctor m_Number;
};

/// provides boundary data for the concentration of the elder problem
/**
 * This class provide boundary data for the elder problem.
 * \todo: This should be temporarily only, until setting of user data can
 * 		  be done efficiently in script/VRL
 */
class ElderConcentrationBoundaryData2d
	: public IBoundaryData<number, 2>
{
	///	Base class type
		typedef IBoundaryData<number, 2> base_type;

	public:
	//	Functor Type
		typedef base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
	///	Constructor
		ElderConcentrationBoundaryData2d()
		{}

	///	virtual destructor
		virtual ~ElderConcentrationBoundaryData2d()	{}

	///	evaluates the data at a given point and time
		bool operator() (number& val, const MathVector<2>& pos, number time = 0.0)
		{
			const number x = pos[0];
			const number y = pos[1];

			if(y == 150){
				if(x > 150 && x < 450){
					val = 1.0;
					return true;
				}
			}
			if(y == 0.0){
				val = 0.0;
				return true;
			}

			val = 0.0;
			return false;
		}
};

/// provides boundary data for the pressure of the elder problem
/**
 * This class provide boundary data for the elder problem.
 * \todo: This should be temporarily only, until setting of user data can
 * 		  be done efficiently in script/VRL
 */
class ElderPressureBoundaryData2d
	: public IBoundaryData<number, 2>
{
	///	Base class type
		typedef IBoundaryData<number, 2> base_type;

	public:
	//	Functor Type
		typedef base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
	///	Constructor
		ElderPressureBoundaryData2d()
		{}

	///	virtual destructor
		virtual ~ElderPressureBoundaryData2d()	{}

	///	evaluates the data at a given point and time
		bool operator() (number& val, const MathVector<2>& pos, number time = 0.0)
		{
			const number x = pos[0];
			const number y = pos[1];

			if(y == 150){
				if(x == 0.0 || x == 600){
					val = 9810 * (150 - y);
					return true;
				}
			}

			val = 0.0;
			return false;
		}
};


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

		virtual void compute(bool compDeriv)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
					value(s, ip) = 1e3 + 0.2e3 * input_value(0, s, ip);
				}

			if(!compDeriv || this->zero_derivative()) return;

			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
					for(size_t fct = 0; fct < num_fct(); ++fct)
						for(size_t dof = 0; dof < num_sh(s, fct); ++dof)
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

		virtual void compute(bool compDeriv)
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
			if(!compDeriv || this->zero_derivative()) return;

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
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
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
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
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
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
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
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
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
						for(size_t sh = 0; sh < num_sh(s, fct); ++sh)
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


template <int dim>
bool RegisterUserData(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	NumberIPData
	{
		typedef IPData<number, dim> T;
		std::stringstream ss; ss << "NumberIPData" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	VectorIPData
	{
		typedef IPData<MathVector<dim>, dim> T;
		std::stringstream ss; ss << "VectorIPData" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	MatrixIPData
	{
		typedef IPData<MathMatrix<dim,dim>, dim> T;
		std::stringstream ss; ss << "MatrixIPData" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	Tensor4IPData
	{
		typedef IPData<MathTensor<4,dim>, dim> T;
		std::stringstream ss; ss << "Tensor4IPData" << dim << "d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str());
	}

//	IUserNumber
	{
		typedef IUserData<number, dim> T;
		typedef IPData<number, dim> TBase;
		std::stringstream ss; ss << "IUserNumber" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str());
	}

//	IUserVector
	{
		typedef IUserData<MathVector<dim>, dim> T;
		typedef IPData<MathVector<dim>, dim> TBase;
		std::stringstream ss; ss << "IUserVector" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str());
	}

//	IUserMatrix
	{
		typedef IUserData<MathMatrix<dim, dim>, dim> T;
		typedef IPData<MathMatrix<dim, dim>, dim> TBase;
		std::stringstream ss; ss << "IUserMatrix" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str());
	}

//	IUserTensor4
	{
		typedef IUserData<MathTensor<4, dim>, dim> T;
		typedef IPData<MathTensor<4, dim>, dim> TBase;
		std::stringstream ss; ss << "IUserTensor4" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str());
	}

//	Base class
	{
		std::stringstream ss; ss << "IBoundaryNumber" << dim << "d";
		reg.add_class_<IBoundaryData<number, dim> >(ss.str().c_str(), grp.c_str());
	}

//	ConstUserNumber
	{
		typedef ConstUserNumber<dim> T;
		std::stringstream ss; ss << "ConstUserNumber" << dim << "d";
		reg.add_class_<T, IPData<number, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set | interactive=false", &T::set, "", "MyNumber || invokeOnChange=true")
			.add_method("print", &T::print);
	}

//	ConstUserVector
	{
		typedef ConstUserVector<dim> T;
		std::stringstream ss; ss << "ConstUserVector" << dim << "d";
		reg.add_class_<T, IUserData<MathVector<dim>, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_all_entries", &T::set_all_entries)
			.add_method("set_entry", &T::set_entry)
			.add_method("print", &T::print);
	}

//	ConstUserMatrix
	{
		typedef ConstUserMatrix<dim> T;
		std::stringstream ss; ss << "ConstUserMatrix" << dim << "d";
		reg.add_class_<T, IUserData<MathMatrix<dim, dim>, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_diag_tensor", &T::set_diag_tensor)
			.add_method("set_all_entries", &T::set_all_entries)
			.add_method("set_entry", &T::set_entry)
			.add_method("print", &T::print);
	}

//	ConstBoundaryNumber
	{
		typedef ConstBoundaryNumber<dim> T;
		std::stringstream ss; ss << "ConstBoundaryNumber" << dim << "d";
		reg.add_class_<T, IBoundaryData<number, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set", &T::set)
			.add_method("print", &T::print);
	}

//	DependentIPDataNumber
	{
		typedef DependentIPData<number, dim> T;
		typedef IPData<number, dim> TBase;
		std::stringstream ss; ss << "DependentIPDataNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DependentIPDataVector
	{
		typedef DependentIPData<MathVector<dim>, dim> T;
		typedef IPData<MathVector<dim>, dim> TBase;
		std::stringstream ss; ss << "DependentIPDataVector" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DependentIPDataMatrix
	{
		typedef DependentIPData<MathMatrix<dim,dim>, dim> T;
		typedef IPData<MathMatrix<dim,dim>, dim> TBase;
		std::stringstream ss; ss << "DependentIPDataMatrix" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DataLinkerNumber
	{
		typedef DataLinker<number, dim> T;
		typedef DependentIPData<number, dim> TBase;
		std::stringstream ss; ss << "DataLinkerNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DataLinkerVector
	{
		typedef DataLinker<MathVector<dim>, dim> T;
		typedef DependentIPData<MathVector<dim>, dim> TBase;
		std::stringstream ss; ss << "DataLinkerVector" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DataLinkerMatrix
	{
		typedef DataLinker<MathMatrix<dim,dim>, dim> T;
		typedef DependentIPData<MathMatrix<dim,dim>, dim> TBase;
		std::stringstream ss; ss << "DataLinkerMatrix" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DataLinkerEqualNumber
	{
		typedef DataLinkerEqualData<number, dim, number> T;
		typedef DataLinker<number, dim> TBase;
		std::stringstream ss; ss << "DataLinkerEqualNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_method("set_input", &T::set_input);
	}

//	DataLinkerEqualMatrix
	{
		typedef DataLinkerEqualData<MathMatrix<dim,dim>, dim, number> T;
		typedef DataLinker<MathMatrix<dim,dim>, dim> TBase;
		std::stringstream ss; ss << "DataLinkerEqualMatrix" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_method("set_input", &T::set_input);
	}

//	IUserFunctionNumber
	{
		typedef IUserFunction<number, dim, number> T;
		typedef DataLinkerEqualData<number, dim, number> TBase;
		std::stringstream ss; ss << "IUserFunctionNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	IUserFunctionMatrixNumber
	{
		typedef IUserFunction<MathMatrix<dim,dim>, dim, number> T;
		typedef DataLinkerEqualData<MathMatrix<dim,dim>, dim, number> TBase;
		std::stringstream ss; ss << "IUserFunctionMatrixNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	ElderDensityLinker
	{
		typedef ElderDensityLinker<dim> T;
		typedef DataLinkerEqualData<number, dim, number> TBase;
		std::stringstream ss; ss << "ElderDensityLinker" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

//	ScaleAddLinkerNumber
	{
		typedef ScaleAddLinker<number, dim, number> T;
		typedef DataLinker<number, dim> TBase;
		std::stringstream ss; ss << "ScaleAddLinkerNumber" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_method("add", &T::add)
			.add_constructor();
	}

//	ScaleAddLinkerVector
	{
		typedef ScaleAddLinker<MathVector<dim>, dim, number> T;
		typedef DataLinker<MathVector<dim>, dim> TBase;
		std::stringstream ss; ss << "ScaleAddLinkerVector" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_method("add", &T::add)
			.add_constructor();
	}

//	ScaleAddLinkerMatrix
	{
		typedef ScaleAddLinker<MathMatrix<dim,dim>, dim, number> T;
		typedef DataLinker<MathMatrix<dim,dim>, dim> TBase;
		std::stringstream ss; ss << "ScaleAddLinkerMatrix" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_method("add", &T::add)
			.add_constructor();
	}

//	DarcyVelocityLinker
	{
		typedef DarcyVelocityLinker<dim> T;
		typedef DataLinker<MathVector<dim>, dim> TBase;
		std::stringstream ss; ss << "DarcyVelocityLinker" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_method("set_density", &T::set_density)
			.add_method("set_gravity", &T::set_gravity)
			.add_method("set_permeability", &T::set_permeability)
			.add_method("set_pressure_gradient", &T::set_pressure_gradient)
			.add_method("set_viscosity", &T::set_viscosity)
			.add_constructor();
	}

	return true;
}

bool RegisterUserData(Registry& reg, const char* parentGroup)
{
//	get group string
	std::string grp = parentGroup; grp.append("/UserData");

#ifdef UG_DIM_1
	RegisterUserData<1>(reg, grp.c_str());
#endif
#ifdef UG_DIM_2
	RegisterUserData<2>(reg, grp.c_str());
#endif
#ifdef UG_DIM_3
	RegisterUserData<3>(reg, grp.c_str());
#endif

#ifdef UG_DIM_2
//	PrintUserNumber2d
//	{
//		typedef PrintUserNumber2d T;
//		std::stringstream ss; ss << "PrintUserNumber2d";
//		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
//			.add_constructor()
//			.add_method("set_user_number|interactive=false", &T::set_user_number, "", "NumberProvider||invokeOnChange=true")
//			.add_method("print", &T::print, "Result", "x#y");
//	}

//	Elder Boundary Data for 2D
	{
		typedef ElderConcentrationBoundaryData2d T;
		typedef IBoundaryData<number, 2> TBase;
		std::stringstream ss; ss << "ElderConcentrationBoundaryData2d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}
	{
		typedef ElderPressureBoundaryData2d T;
		typedef IBoundaryData<number, 2> TBase;
		std::stringstream ss; ss << "ElderPressureBoundaryData2d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
			.add_constructor();
	}

#endif

	return true;
}

} // end namespace
} // end namepace

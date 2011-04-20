

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
	: public DataLinker<number, dim, number>
{
	///	Base class type
		typedef DataLinker<number, dim, number> base_type;

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


template <int dim>
bool RegisterUserData(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	Base class
	{
		std::stringstream ss; ss << "NumberIPData" << dim << "d";
		reg.add_class_<IPData<number, dim> >(ss.str().c_str(), grp.c_str());
	}

//	Base class
	{
		std::stringstream ss; ss << "VectorIPData" << dim << "d";
		reg.add_class_<IPData<MathVector<dim>, dim> >(ss.str().c_str(), grp.c_str());
	}

//	Base class
	{
		std::stringstream ss; ss << "MatrixIPData" << dim << "d";
		reg.add_class_<IPData<MathMatrix<dim,dim>, dim> >(ss.str().c_str(), grp.c_str());
	}

//	Base class
	{
		std::stringstream ss; ss << "IUserNumber" << dim << "d";
		reg.add_class_<IUserData<number, dim>, IPData<number, dim> >(ss.str().c_str(), grp.c_str());
	}

//	Base class
	{
		std::stringstream ss; ss << "IUserVector" << dim << "d";
		reg.add_class_<IUserData<MathVector<dim>, dim>, IPData<MathVector<dim>, dim> >(ss.str().c_str(), grp.c_str());
	}

//	Base class
	{
		std::stringstream ss; ss << "IUserMatrix" << dim << "d";
		reg.add_class_<IUserData<MathMatrix<dim, dim>, dim>, IPData<MathMatrix<dim, dim>, dim> >(ss.str().c_str(), grp.c_str());
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

//	DependentIPDataMatrix
	{
		typedef DependentIPData<MathMatrix<dim,dim>, dim> T;
		typedef IPData<MathMatrix<dim,dim>, dim> TBase;
		std::stringstream ss; ss << "DependentIPDataMatrix" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	DataLinkerNumber
	{
		typedef DataLinker<number, dim, number> T;
		typedef DependentIPData<number, dim> TBase;
		std::stringstream ss; ss << "DataLinkerNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_method("set_input", &T::set_input);
	}

//	DataLinkerMatrixNumber
	{
		typedef DataLinker<MathMatrix<dim,dim>, dim, number> T;
		typedef DependentIPData<MathMatrix<dim,dim>, dim> TBase;
		std::stringstream ss; ss << "DataLinkerMatrixNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str())
			.add_method("set_input", &T::set_input);
	}

//	IUserFunctionNumber
	{
		typedef IUserFunction<number, dim, number> T;
		typedef DataLinker<number, dim, number> TBase;
		std::stringstream ss; ss << "IUserFunctionNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	IUserFunctionMatrixNumber
	{
		typedef IUserFunction<MathMatrix<dim,dim>, dim, number> T;
		typedef DataLinker<MathMatrix<dim,dim>, dim, number> TBase;
		std::stringstream ss; ss << "IUserFunctionMatrixNumber" << dim << "d";
		reg.add_class_<T, TBase >(ss.str().c_str(), grp.c_str());
	}

//	ElderDensityLinker
	{
		typedef ElderDensityLinker<dim> T;
		typedef DataLinker<number, dim, number> TBase;
		std::stringstream ss; ss << "ElderDensityLinker" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), grp.c_str())
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

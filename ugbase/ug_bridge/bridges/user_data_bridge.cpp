

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"
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
		reg.add_class_<IBoundaryNumberProvider<dim> >(ss.str().c_str(), grp.c_str());
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
		reg.add_class_<T, IBoundaryNumberProvider<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set", &T::set)
			.add_method("print", &T::print);
	}

	return true;
}

bool RegisterUserData(Registry& reg, const char* parentGroup)
{
//	get group string
	std::string grp = parentGroup; grp.append("/UserData");

	RegisterUserData<1>(reg, grp.c_str());
	RegisterUserData<2>(reg, grp.c_str());
	RegisterUserData<3>(reg, grp.c_str());

//	PrintUserNumber2d
	{
		typedef PrintUserNumber2d T;
		std::stringstream ss; ss << "PrintUserNumber2d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_user_number|interactive=false", &T::set_user_number, "", "NumberProvider||invokeOnChange=true")
			.add_method("print", &T::print, "Result", "x#y");
	}

	return true;
}

} // end namespace
} // end namepace

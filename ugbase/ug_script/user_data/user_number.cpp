

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/spacial_discretization/user_data.h"
#include "ug_script/ug_script.h"
#include <iostream>
#include <sstream>

namespace ug
{
namespace bridge
{

template <int dim>
class LuaUserNumber : public IUserNumberProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserNumberProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
		LuaUserNumber()
		{
			m_L = ug::script::GetDefaultLuaState();
		}

		void set_lua_callback(const char* luaCallback)
		{
		//	Todo: Work on function pointer directly
			m_callbackName = luaCallback;
		}

		void operator() (number& c, const MathVector<dim>& x, number time = 0.0) const
		{
			lua_getglobal(m_L, m_callbackName);
			for(size_t i = 0; i < dim; ++i)
				lua_pushnumber(m_L, x[i]);
			lua_pushnumber(m_L, time);

			if(lua_pcall(m_L, dim + 1, 1, 0) != 0)
			{
				UG_LOG("error running lua number callback '" << m_callbackName << "': "
								<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

			c = luaL_checknumber(m_L, -1);
			lua_pop(m_L, 1);
		}

	protected:
		const char* m_callbackName;
		lua_State*	m_L;
};

class PrintUserNumber2d
{
	protected:
		typedef IUserNumberProvider<2>::functor_type NumberFunctor;

	public:
		void set_user_number(IUserNumberProvider<2>& user)
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
void RegisterUserNumber(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	Base class
	{
		std::stringstream ss; ss << "IUserNumberProvider" << dim << "d";
		reg.add_class_<IUserNumberProvider<dim> >(ss.str().c_str(), grp.c_str());
	}

//	ConstUserNumber
	{
		typedef ConstUserNumber<dim> T;
		std::stringstream ss; ss << "ConstUserNumber" << dim << "d";
		reg.add_class_<T, IUserNumberProvider<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set | interactive=false", &T::set, "", "MyNumber || invokeOnChange=true")
			.add_method("print", &T::print);
	}

//	LuaUserNumber
	{
		typedef LuaUserNumber<dim> T;
		std::stringstream ss; ss << "LuaUserNumber" << dim << "d";
		reg.add_class_<T, IUserNumberProvider<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}
}

void RegisterUserNumber(Registry& reg, const char* parentGroup)
{
	RegisterUserNumber<1>(reg, parentGroup);
	RegisterUserNumber<2>(reg, parentGroup);
	RegisterUserNumber<3>(reg, parentGroup);

//	PrintUserNumber2d
	{
		std::string grp = std::string(parentGroup);

		typedef PrintUserNumber2d T;
		std::stringstream ss; ss << "PrintUserNumber2d";
		reg.add_class_<T>(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_user_number|interactive=false", &T::set_user_number, "", "NumberProvider||invokeOnChange=true")
			.add_method("print", &T::print, "Result", "x#y");
	}

}

} // end namespace
} // end namepace

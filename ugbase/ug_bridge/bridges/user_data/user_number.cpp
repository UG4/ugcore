

#include "../../ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/lib_discretization.h"
#include "../../../ug_script/ug_script.h"

namespace ug
{
namespace bridge
{

template <int dim>
class ConstUserNumber
{
	public:
		ConstUserNumber() {m_Number = 0.0;}

		void set(number val)
		{
			m_Number = val;
		}

		void print() const
		{
			UG_LOG("ConstUserNumber:" << m_Number << "\n");
		}

		void operator() (number& c, const MathVector<dim>& x, number time = 0.0)
		{
			c = m_Number;
		}

	protected:
		number m_Number;
};

template <int dim>
class LuaUserNumber
{
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

		void operator() (number& c, const MathVector<dim>& x, number time = 0.0)
		{
			lua_getglobal(m_L, m_callbackName);
			for(size_t i = 0; i < dim; ++i)
				lua_pushnumber(m_L, x[i]);

			if(lua_pcall(m_L, dim, 1, 0) != 0)
			{
				UG_LOG("error running diffusion callback " << m_callbackName << ": "
								<< lua_tostring(m_L, -1));
				throw(int(0));
			}

			c = luaL_checknumber(m_L, -1);
			lua_pop(m_L, 1);
		}

	protected:
		const char* m_callbackName;
		lua_State*	m_L;
};


template <int dim, typename TUserNumber>
class UserNumberProvider : public IUserNumberProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserNumberProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return m_UserNumber;}

	//	set user Number
		void set_user_number(const TUserNumber& userNumber) {m_UserNumber = userNumber;}

	protected:
		TUserNumber	m_UserNumber;
};


void RegisterUserNumber(Registry& reg)
{
	static const int dim = 2;

//	UserNumberProvider
	{
	//	Base class
		reg.add_class_<IUserNumberProvider<dim> >("IUserNumberProvider2d");

	//	Const Number Provider
		{
			typedef UserNumberProvider<dim, ConstUserNumber<dim> > T;
			reg.add_class_<T, IUserNumberProvider<dim> >("ConstUserNumberProvider2d")
				.add_constructor()
				.add_method("set_user_number", &T::set_user_number);
		}

	//	Lua Number
		{
			typedef UserNumberProvider<dim, LuaUserNumber<dim> > T;
			reg.add_class_<T, IUserNumberProvider<dim> >("LuaUserNumberProvider2d")
				.add_constructor()
				.add_method("set_user_number", &T::set_user_number);
		}
	}

//	Functor
	{
	//	ConstUserNumber
		{
			typedef ConstUserNumber<dim> T;
			reg.add_class_<T>("ConstUserNumber2d")
				.add_constructor()
				.add_method("set", &T::set)
				.add_method("print", &T::print);
		}

	//	LuaUserNumber
		{
			typedef LuaUserNumber<dim> T;
			reg.add_class_<T>("LuaUserNumber2d")
				.add_constructor()
				.add_method("set_lua_callback", &T::set_lua_callback);
		}
	}
}

} // end namespace
} // end namepace

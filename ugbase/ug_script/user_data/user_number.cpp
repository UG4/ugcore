

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/spatial_discretization/user_data.h"
#include "ug_script/ug_script.h"
#include <iostream>
#include <sstream>

namespace ug
{
namespace bridge
{

template <int dim>
class LuaUserNumber : public IUserNumber<dim>
{
	public:
	//	Functor Type
		typedef typename IUserNumber<dim>::functor_type functor_type;

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

template <int dim>
void RegisterLuaUserNumber(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	LuaUserNumber
	{
		typedef LuaUserNumber<dim> T;
		std::stringstream ss; ss << "LuaUserNumber" << dim << "d";
		reg.add_class_<T, IUserNumber<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}
}

void RegisterLuaUserNumber(Registry& reg, const char* parentGroup)
{
	RegisterLuaUserNumber<1>(reg, parentGroup);
	RegisterLuaUserNumber<2>(reg, parentGroup);
	RegisterLuaUserNumber<3>(reg, parentGroup);
}

} // end namespace
} // end namepace

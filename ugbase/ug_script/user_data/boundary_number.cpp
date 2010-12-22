

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

template <int dim>
class LuaBoundaryNumber : public IBoundaryNumberProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IBoundaryNumberProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
		LuaBoundaryNumber()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

		virtual ~LuaBoundaryNumber()	{}

		void set_lua_callback(const char* luaCallback)
		{
			m_callbackName = luaCallback;
		//	store the callback function in the registry and obtain a reference.
			lua_getglobal(m_L, m_callbackName);
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

		bool operator() (number& c, const MathVector<dim>& x, number time = 0.0)
		{
		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

			for(size_t i = 0; i < dim; ++i)
				lua_pushnumber(m_L, x[i]);
			lua_pushnumber(m_L, time);

			if(lua_pcall(m_L, dim + 1, 2, 0) != 0)
			{
				UG_LOG("error running diffusion callback " << m_callbackName << ": "
								<< lua_tostring(m_L, -1));
				throw(int(0));
			}

			c = luaL_checknumber(m_L, -1);
			bool res = lua_toboolean(m_L, -2);
			lua_pop(m_L, 2);
			return res;
		}

	protected:
		const char* m_callbackName;
		int m_callbackRef;
		lua_State*	m_L;
};

template <int dim>
void RegisterLuaBoundaryNumber(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	LuaBoundaryNumber
	{
		typedef LuaBoundaryNumber<dim> T;
		std::stringstream ss; ss << "LuaBoundaryNumber" << dim << "d";
		reg.add_class_<T, IBoundaryNumberProvider<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}
}

void RegisterLuaBoundaryNumber(Registry& reg, const char* parentGroup)
{
	RegisterLuaBoundaryNumber<1>(reg, parentGroup);
	RegisterLuaBoundaryNumber<2>(reg, parentGroup);
	RegisterLuaBoundaryNumber<3>(reg, parentGroup);
}

} // end namespace
} // end namepace

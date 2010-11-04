

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "lib_discretization/spacial_discretization/user_data.h"
#include "ug_script/ug_script.h"
#include <sstream>
#include <string>

namespace ug
{
namespace bridge
{

template <int dim>
class LuaUserVector : public IUserVectorProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserVectorProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
		LuaUserVector()
		{
			m_L = ug::script::GetDefaultLuaState();
		}

		void set_lua_callback(const char* luaCallback)
		{
		//	Todo: Work on function pointer directly
			m_callbackName = luaCallback;
		}

		void operator() (MathVector<dim>& v, const MathVector<dim>& x, number time = 0.0)
		{
			lua_getglobal(m_L, m_callbackName);
			for(size_t i = 0; i < dim; ++i)
				lua_pushnumber(m_L, x[i]);
			lua_pushnumber(m_L, time);

			if(lua_pcall(m_L, dim + 1, dim, 0) != 0)
			{
				UG_LOG("error running lua vector callback '" << m_callbackName << "': "
								<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

			int counter = -1;
			for(size_t i = 0; i < dim; ++i){
					v[i] = luaL_checknumber(m_L, counter--);
			}

			lua_pop(m_L, dim);
		}

	protected:
		const char* m_callbackName;
		lua_State*	m_L;
};


template <int dim>
void RegisterLuaUserVector(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

	//	LuaUserVector
	{
		typedef LuaUserVector<dim> T;
		std::stringstream ss; ss << "LuaUserVector" << dim << "d";
		reg.add_class_<T, IUserVectorProvider<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}
}

void RegisterLuaUserVector(Registry& reg, const char* parentGroup)
{
	RegisterLuaUserVector<1>(reg, parentGroup);
	RegisterLuaUserVector<2>(reg, parentGroup);
	RegisterLuaUserVector<3>(reg, parentGroup);
}

} // end namespace
} // end namepace

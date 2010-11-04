

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
class LuaUserMatrix : public IUserMatrixProvider<dim>
{
	public:
	//	Functor Type
		typedef typename IUserMatrixProvider<dim>::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
		LuaUserMatrix()
		{
			m_L = ug::script::GetDefaultLuaState();
		}

		void set_lua_callback(const char* luaCallback)
		{
		//	Todo: Work on function pointer directly
			m_callbackName = luaCallback;
		}

		void operator() (MathMatrix<dim, dim>& D, const MathVector<dim>& x, number time = 0.0)
		{
			lua_getglobal(m_L, m_callbackName);
			for(size_t i = 0; i < dim; ++i)
				lua_pushnumber(m_L, x[i]);
			lua_pushnumber(m_L, time);

			size_t matSize = dim*dim;

			if(lua_pcall(m_L, dim + 1, matSize, 0) != 0)
			{
				UG_LOG("error running lua matrix callback '" << m_callbackName << "': "
								<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

			int counter = -1;
			for(size_t i = 0; i < dim; ++i){
				for(size_t j = 0; j < dim; ++j){
					D[i][j] = luaL_checknumber(m_L, counter--);
				}
			}

			lua_pop(m_L, matSize);
		}

	protected:
		const char* m_callbackName;
		lua_State*	m_L;
};

template <int dim>
void RegisterLuaUserMatrix(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	LuaUserMatrix
	{
		typedef LuaUserMatrix<dim> T;
		std::stringstream ss; ss << "LuaUserMatrix" << dim << "d";
		reg.add_class_<T, IUserMatrixProvider<dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}
}

void RegisterLuaUserMatrix(Registry& reg, const char* parentGroup)
{
	RegisterLuaUserMatrix<1>(reg, parentGroup);
	RegisterLuaUserMatrix<2>(reg, parentGroup);
	RegisterLuaUserMatrix<3>(reg, parentGroup);
}


} // end namespace
} // end namepace

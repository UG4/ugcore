

#include "ug_bridge/ug_bridge.h"
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_discretization/spatial_discretization/ip_data/user_data.h"
#include "ug_script/ug_script.h"
#include <iostream>
#include <sstream>
#include "./user_data.h"

namespace ug
{
namespace bridge
{

/// Lua Traits to push/pop on lua stack
template <typename TData>
struct lua_traits;

template <>
struct lua_traits<number>
{
	static void push(lua_State*	L, const number& c)
	{
		lua_pushnumber(L, c);
	}

	static void read(lua_State* L, number& c)
	{
		c = luaL_checknumber(L, -1);
	}

	static const int size = 1;
};

template <std::size_t dim>
struct lua_traits< ug::MathVector<dim> >
{
	static void push(lua_State*	L, const MathVector<dim>& x)
	{
		for(size_t i = 0; i < dim; ++i)
			lua_pushnumber(L, x[i]);
	}

	static void read(lua_State* L, MathVector<dim>& x)
	{
		int counter = -1;
		for(size_t i = 0; i < dim; ++i){
				x[dim-1-i] = luaL_checknumber(L, counter--);
		}
	}

	static const int size = dim;
};


template <std::size_t dim>
struct lua_traits< MathMatrix<dim, dim> >
{
	static void push(lua_State*	L, const MathMatrix<dim, dim>& D)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				lua_pushnumber(L, D[i][j]);
			}
		}

	}

	static void read(lua_State* L, MathMatrix<dim, dim>& D)
	{
		int counter = -1;
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				D[dim-1-j][dim-1-i] = luaL_checknumber(L, counter--);
			}
		}
	}

	static const int size = dim*dim;
};

////////////////////////////////
// Generic LuaUserData
////////////////////////////////

template <typename TData, int dim>
class LuaUserData
	: public IUserData<TData, dim>
{
	///	Base class type
		typedef IUserData<TData, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	//	Functor Type
		typedef typename base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
		LuaUserData()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

		virtual ~LuaUserData()	{}

		void set_lua_callback(const char* luaCallback)
		{
			m_callbackName = luaCallback;
		//	store the callback function in the registry and obtain a reference.
			lua_getglobal(m_L, m_callbackName);
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

		void operator() (TData& D, const MathVector<dim>& x, number time = 0.0)
		{
		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

		//  push space coordinates on stack
			lua_traits<MathVector<dim> >::push(m_L, x);

		//	push time on stack
			lua_traits<number>::push(m_L, time);

		//	compute total args size
			size_t argSize = lua_traits<MathVector<dim> >::size
								+  lua_traits<number>::size;

		//	compute total return size
			size_t retSize = lua_traits<TData>::size;

		//	call lua function
			if(lua_pcall(m_L, argSize, retSize, 0) != 0)
			{
				UG_LOG("error running lua matrix callback '" << m_callbackName << "': "
								<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

		//	read return value
			lua_traits<TData>::read(m_L, D);

		//	pop values
			lua_pop(m_L, retSize);
		}

	///	implement as a IPData
		virtual void compute(bool computeDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
				{
					this->operator()(	value(s,i),
										ip(s, i),
										time());
				}
		}

	protected:
		const char* m_callbackName;
		int m_callbackRef;
		lua_State*	m_L;
};

////////////////////////////////
// Generic LuaUserFunction
////////////////////////////////

/// this class maps a scalar value an output scalar value using a lua callback
template <typename TData>
class LuaUserFunction
	: public IUserFunction<TData>
{
	private:
	//	type of base class
		typedef IUserFunction<TData> base_type;

	public:
		LuaUserFunction()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

		void set_lua_callback(const char* luaCallback)
		{
			m_callbackName = luaCallback;
		//	store the callback function in the registry and obtain a reference.
			lua_getglobal(m_L, m_callbackName);
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

		void operator() (TData& out, int numArgs, ...)
		{

		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

			va_list ap;
			va_start(ap, numArgs);

			for(int i = 0; i < numArgs; ++i)
			{
				number val = va_arg(ap, number);
				lua_pushnumber(m_L, val);
				UG_LOG("Push value i=" << i << ": " << val<<"\n");
			}

			va_end(ap);


			if(lua_pcall(m_L, numArgs, 1, 0) != 0)
			{
				UG_LOG("error running lua number callback '" << m_callbackName << "': "
								<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

			out = luaL_checknumber(m_L, -1);
			lua_pop(m_L, 1);
		}

	protected:
		const char* m_callbackName;
		int m_callbackRef;
		lua_State*	m_L;
};

LuaUserNumberNumberFunction::LuaUserNumberNumberFunction()
{
	m_L = ug::script::GetDefaultLuaState();
	m_callbackRef = LUA_NOREF;
}

void LuaUserNumberNumberFunction::set_lua_callback(const char* luaCallback)
{
	m_callbackName = luaCallback;
//	store the callback function in the registry and obtain a reference.
	lua_getglobal(m_L, m_callbackName);
	m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
}

number LuaUserNumberNumberFunction::operator() (int numArgs, ...) const
{

//	push the callback function on the stack
	lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

	va_list ap;
	va_start(ap, numArgs);

	for(int i = 0; i < numArgs; ++i)
	{
		number val = va_arg(ap, number);
		lua_pushnumber(m_L, val);
		UG_LOG("Push value i=" << i << ": " << val<<"\n");
	}

	va_end(ap);


	if(lua_pcall(m_L, numArgs, 1, 0) != 0)
	{
		UG_LOG("error running lua number callback '" << m_callbackName << "': "
						<< lua_tostring(m_L, -1) << "\n");
		throw(int(0));
	}

	number c = luaL_checknumber(m_L, -1);
	lua_pop(m_L, 1);

	return c;
}

template <int dim>
void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	LuaUserNumber
	{
		typedef LuaUserData<number, dim> T;
		std::stringstream ss; ss << "LuaUserNumber" << dim << "d";
		reg.add_class_<T, IUserData<number, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaUserVector
	{
		typedef LuaUserData<MathVector<dim>, dim> T;
		std::stringstream ss; ss << "LuaUserVector" << dim << "d";
		reg.add_class_<T, IUserData<MathVector<dim>, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaUserMatrix
	{
		typedef LuaUserData<MathMatrix<dim,dim>, dim> T;
		std::stringstream ss; ss << "LuaUserMatrix" << dim << "d";
		reg.add_class_<T, IUserData<MathMatrix<dim, dim>, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

}

void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{

//	LuaUserFunction
	{
		typedef LuaUserFunction<number> T;
		reg.add_class_<T>("LuaUserFunctionNumber", parentGroup)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaUserNumberNumberFunction
	{
		typedef LuaUserNumberNumberFunction T;
		std::stringstream ss; ss << "LuaUserNumberNumberFunction";
		reg.add_class_<T>(ss.str().c_str(), parentGroup)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

	RegisterLuaUserData<1>(reg, parentGroup);
	RegisterLuaUserData<2>(reg, parentGroup);
	RegisterLuaUserData<3>(reg, parentGroup);
}

} // end namespace
} // end namepace

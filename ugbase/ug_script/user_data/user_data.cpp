

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
class LuaUserNumber
	: public IUserData<number, dim>
{
	/// Base class type
		typedef IUserData<number, dim> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	//	Functor Type
		typedef typename base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return boost::ref(*this);}

	public:
		LuaUserNumber()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

		virtual ~LuaUserNumber()	{}

		void set_lua_callback(const char* luaCallback)
		{
			m_callbackName = luaCallback;
		//	store the callback function in the registry and obtain a reference.
			lua_getglobal(m_L, m_callbackName);
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

		void operator() (number& c, const MathVector<dim>& x, number time = 0.0) const
		{
		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

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

template <int dim>
class LuaUserVector
	: public IUserData<MathVector<dim>, dim>
{
/// Base class type
	typedef IUserData<MathVector<dim>, dim> base_type;

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
		LuaUserVector()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

		virtual ~LuaUserVector()	{}

		void set_lua_callback(const char* luaCallback)
		{
			m_callbackName = luaCallback;
		//	store the callback function in the registry and obtain a reference.
			lua_getglobal(m_L, m_callbackName);
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

		void operator() (MathVector<dim>& v, const MathVector<dim>& x, number time = 0.0)
		{
		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

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
					v[dim-1-i] = luaL_checknumber(m_L, counter--);
			}

			lua_pop(m_L, dim);
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


template <int dim>
class LuaUserMatrix
	: public IUserData<MathMatrix<dim, dim>, dim>
{
	///	Base class type
		typedef IUserData<MathMatrix<dim, dim>, dim> base_type;

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
		LuaUserMatrix()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

		virtual ~LuaUserMatrix()	{}

		void set_lua_callback(const char* luaCallback)
		{
			m_callbackName = luaCallback;
		//	store the callback function in the registry and obtain a reference.
			lua_getglobal(m_L, m_callbackName);
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

		void operator() (MathMatrix<dim, dim>& D, const MathVector<dim>& x, number time = 0.0)
		{
		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

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
					D[dim-1-j][dim-1-i] = luaL_checknumber(m_L, counter--);
				}
			}

			lua_pop(m_L, matSize);
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

template <int dim>
void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{
	std::string grp = std::string(parentGroup);

//	LuaUserNumber
	{
		typedef LuaUserNumber<dim> T;
		std::stringstream ss; ss << "LuaUserNumber" << dim << "d";
		reg.add_class_<T, IUserData<number, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaUserVector
	{
		typedef LuaUserVector<dim> T;
		std::stringstream ss; ss << "LuaUserVector" << dim << "d";
		reg.add_class_<T, IUserData<MathVector<dim>, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaUserMatrix
	{
		typedef LuaUserMatrix<dim> T;
		std::stringstream ss; ss << "LuaUserMatrix" << dim << "d";
		reg.add_class_<T, IUserData<MathMatrix<dim, dim>, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

}

void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{
	RegisterLuaUserData<1>(reg, parentGroup);
	RegisterLuaUserData<2>(reg, parentGroup);
	RegisterLuaUserData<3>(reg, parentGroup);
}

} // end namespace
} // end namepace

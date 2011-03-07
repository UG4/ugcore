

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

/// provides data specified in the lua script
/**
 * This class implements the IPData interface to provide data at arbitrary
 * points. A Lua callback is used to evaluate the data.
 */
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
	///	Constructor
		LuaUserData()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

	///	virtual destructor
		virtual ~LuaUserData()	{}

	///	sets the Lua function used to compute the data
	/**
	 * This function sets the lua callback. The name of the function is
	 * passed as a string. Make sure, that the function name is defined
	 * when executing the script.
	 */
		void set_lua_callback(const char* luaCallback)
		{
		//	store name (string) of callback
			m_callbackName = luaCallback;

		//	obtain a reference
			lua_getglobal(m_L, m_callbackName);

		//	store reference to lua function
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

	///	evaluates the data at a given point and time
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
				UG_LOG("ERROR in 'LuaUserData::operator(...)': Error while "
						"running lua matrix callback '" << m_callbackName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n");
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
	///	callback name as string
		const char* m_callbackName;

	///	reference to lua function
		int m_callbackRef;

	///	lua state
		lua_State*	m_L;
};

////////////////////////////////
// Generic LuaBoundaryData
////////////////////////////////

/// provides boundary data specified in the lua script
/**
 * This class provide boundary data at arbitrary
 * points. A Lua callback is used to evaluate the data.
 */
template <typename TData, int dim>
class LuaBoundaryData
	: public IBoundaryData<TData, dim>
{
	///	Base class type
		typedef IBoundaryData<TData, dim> base_type;

	public:
	//	Functor Type
		typedef typename base_type::functor_type functor_type;

	//	return functor
		virtual functor_type get_functor() const {return *this;}

	public:
	///	Constructor
		LuaBoundaryData()
		{
			m_L = ug::script::GetDefaultLuaState();
			m_callbackRef = LUA_NOREF;
		}

	///	virtual destructor
		virtual ~LuaBoundaryData()	{}

	///	sets the Lua function used to compute the data
	/**
	 * This function sets the lua callback. The name of the function is
	 * passed as a string. Make sure, that the function name is defined
	 * when executing the script.
	 */
		void set_lua_callback(const char* luaCallback)
		{
		//	store name (string) of callback
			m_callbackName = luaCallback;

		//	obtain a reference
			lua_getglobal(m_L, m_callbackName);

		//	store reference to lua function
			m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
		}

	///	evaluates the data at a given point and time
		bool operator() (TData& D, const MathVector<dim>& x, number time = 0.0)
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

		//	compute total return size (+1 for boolian)
			size_t retSize = lua_traits<TData>::size + 1;

		//	call lua function
			if(lua_pcall(m_L, argSize, retSize, 0) != 0)
			{
				UG_LOG("ERROR in 'LuaUserData::operator(...)': Error while "
						"running lua matrix callback '" << m_callbackName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

		//	read return value
			lua_traits<TData>::read(m_L, D);

		//	read bool flag
			bool res = lua_toboolean(m_L, -retSize);

		//	pop values
			lua_pop(m_L, retSize);

		//	forward flag
			return res;
		}

	protected:
	///	callback name as string
		const char* m_callbackName;

	///	reference to lua function
		int m_callbackRef;

	///	lua state
		lua_State*	m_L;
};

////////////////////////////////
// Generic LuaUserFunction
////////////////////////////////

/// maps several data values to an output data value using a lua callback
/**
 * This class provides the evaluation of a user function, that is specified
 * in the script. Several data (of the same c++-type) can be used as input,
 * a data (of the same type) is returned.
 */
template <typename TData, int dim, typename TDataIn>
class LuaUserFunction
	: public IUserFunction<TData, dim, TDataIn>
{
	private:
	//	type of base class
		typedef IUserFunction<TData, dim, TDataIn> base_type;

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
		using base_type::num_input;
		using base_type::input_value;
		using base_type::input_deriv;

	public:
	///	constructor
		LuaUserFunction() : m_numArgs(0)
		{
			m_L = ug::script::GetDefaultLuaState();
			m_cbValueRef = LUA_NOREF;
			m_cbDerivRef.clear();
			m_cbDerivName.clear();
		}

	///	sets the Lua function used to compute the data
	/**
	 * This function sets the lua callback. The name of the function is
	 * passed as a string. Make sure, that the function name is defined
	 * when executing the script.
	 */
		void set_lua_value_callback(const char* luaCallback, size_t numArgs)
		{
		//	store name (string) of callback
			m_cbValueName = luaCallback;

		//	obtain a reference
			lua_getglobal(m_L, m_cbValueName);

		//	store reference to lua function
			m_cbValueRef = luaL_ref(m_L, LUA_REGISTRYINDEX);

		//	remember number of arguments to be used
			m_numArgs = numArgs;
			m_cbDerivName.resize(numArgs);
			m_cbDerivRef.resize(numArgs, LUA_NOREF);

		//	set num inputs for linker
			set_num_input(numArgs);
		}

	///	sets the Lua function used to compute the derivative
		bool set_lua_deriv_callback(size_t arg, const char* luaCallback)
		{
		//	check number of arg
			if(arg >= m_numArgs)
			{
				UG_LOG("ERROR in set_lua_deriv_callback: Trying to set a derivative"
						" for an argument, that does not exist.\n");
				return false;
			}

		//	store name (string) of callback
			m_cbDerivName[arg] = luaCallback;

		//	obtain a reference
			lua_getglobal(m_L, m_cbDerivName[arg]);

		//	store reference to lua function
			m_cbDerivRef[arg] = luaL_ref(m_L, LUA_REGISTRYINDEX);

			return true;
		}

	///	evaluates the data
		virtual void operator() (TData& out, int numArgs, ...)
		{
			UG_ASSERT(numArgs == (int)m_numArgs, "Number of arguments mismatched.");

		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_cbValueRef);

		//	get list of arguments
			va_list ap;
			va_start(ap, numArgs);

		//	read all arguments and push them to the lua stack
			for(int i = 0; i < numArgs; ++i)
			{
			//	cast data
				TDataIn val = va_arg(ap, TDataIn);

			//	push data to lua stack
				lua_traits<TDataIn>::push(m_L, val);
			}

		//	end read in of parameters
			va_end(ap);

		//	compute total args size
			size_t argSize = lua_traits<TDataIn>::size * numArgs;

		//	compute total return size
			size_t retSize = lua_traits<TData>::size;

		//	call lua function
			if(lua_pcall(m_L, argSize, retSize, 0) != 0)
			{
				UG_LOG("ERROR in 'LuaUserFunction::operator(...)': Error while "
						"running lua matrix callback '" << m_cbValueName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

		//	read return value
			lua_traits<TData>::read(m_L, out);

		//	pop values
			lua_pop(m_L, retSize);
		}

	///	computes the value
		virtual void compute(bool compDeriv)
		{
		//	vector of data for all inputs
			std::vector<TDataIn> vDataIn(num_input());

			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
				//	gather all input data for this ip
					for(size_t c = 0; c < vDataIn.size(); ++c)
						vDataIn[c] = input_value(c, s, ip);

				//	evaluate data at ip
					eval_value(value(s,ip), vDataIn);
				}

		//	check if derivative is required
			if(!compDeriv || this->zero_derivative()) return;

			for(size_t s = 0; s < num_series(); ++s)
				for(size_t ip = 0; ip < num_ip(s); ++ip)
				{
					for(size_t c = 0; c < vDataIn.size(); ++c)
					{
					//	gather all input data for this ip
						vDataIn[c] = input_value(c, s, ip);

					//	data of derivative w.r.t. one component at ip-values
						TData derivVal;

					//	evaluate data at ip
						eval_deriv(derivVal, vDataIn, c);

					//	loop functions
						for(size_t fct = 0; fct < num_fct(); ++fct)
							for(size_t dof = 0; dof < num_sh(s, fct); ++dof)
							{
							//	todo: Is this product always defined?
								deriv(s, ip, fct, dof) = derivVal * input_deriv(0, s, ip, fct, dof);
							}
					}
				}
		}

	///	evaluates the data
		void eval_value(TData& out, const std::vector<TDataIn>& dataIn)
		{
			UG_ASSERT(dataIn.size() == m_numArgs, "Number of arguments mismatched.");

		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_cbValueRef);

		//	read all arguments and push them to the lua stack
			for(size_t i = 0; i < dataIn.size(); ++i)
			{
			//	push data to lua stack
				lua_traits<TDataIn>::push(m_L, dataIn[i]);
			}

		//	compute total args size
			size_t argSize = lua_traits<TDataIn>::size * dataIn.size();

		//	compute total return size
			size_t retSize = lua_traits<TData>::size;

		//	call lua function
			if(lua_pcall(m_L, argSize, retSize, 0) != 0)
			{
				UG_LOG("ERROR in 'LuaUserFunction::operator(...)': Error while "
						"running lua matrix callback '" << m_cbValueName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

		//	read return value
			lua_traits<TData>::read(m_L, out);

		//	pop values
			lua_pop(m_L, retSize);
		}


	///	evaluates the data
		void eval_deriv(TData& out, const std::vector<TDataIn>& dataIn, size_t arg)
		{
			UG_ASSERT(dataIn.size() == m_numArgs, "Number of arguments mismatched.");
			UG_ASSERT(arg < m_numArgs, "Argument does not exist.");

		//	push the callback function on the stack
			lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_cbDerivRef[arg]);

		//	read all arguments and push them to the lua stack
			for(size_t i = 0; i < dataIn.size(); ++i)
			{
			//	push data to lua stack
				lua_traits<TDataIn>::push(m_L, dataIn[i]);
			}

		//	compute total args size
			size_t argSize = lua_traits<TDataIn>::size * dataIn.size();

		//	compute total return size
			size_t retSize = lua_traits<TData>::size;

		//	call lua function
			if(lua_pcall(m_L, argSize, retSize, 0) != 0)
			{
				UG_LOG("ERROR in 'LuaUserFunction::operator(...)': Error while "
						"running lua matrix callback '" << m_cbDerivName[arg] << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n");
				throw(int(0));
			}

		//	read return value
			lua_traits<TData>::read(m_L, out);

		//	pop values
			lua_pop(m_L, retSize);
		}

	protected:
	///	callback name as string
		const char* m_cbValueName;
		std::vector<const char*> m_cbDerivName;

	///	reference to lua function
		int m_cbValueRef;
		std::vector<int> m_cbDerivRef;

	///	lua state
		lua_State*	m_L;

	///	number of arguments to use
		size_t m_numArgs;
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

//	LuaBoundaryNumber
	{
		typedef LuaBoundaryData<number, dim> T;
		std::stringstream ss; ss << "LuaBoundaryNumber" << dim << "d";
		reg.add_class_<T, IBoundaryData<number, dim> >(ss.str().c_str(), grp.c_str())
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaUserFunction
	{
		typedef LuaUserFunction<number, dim, number> T;
		typedef IUserFunction<number, dim, number> TBase;
		std::stringstream ss; ss << "LuaUserFunctionNumber" << dim << "d";
		reg.add_class_<T, TBase>(ss.str().c_str(), parentGroup)
			.add_constructor()
			.add_method("set_lua_value_callback", &T::set_lua_value_callback)
			.add_method("set_lua_deriv_callback", &T::set_lua_deriv_callback);
	}
}

void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{

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

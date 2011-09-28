// author: Andreas Vogel

#include <iostream>
#include <sstream>

#include "bridge/ug_bridge.h"
#include "common/common.h"
#include "common/math/ugmath.h"
#include "lua_user_data.h"
#include "lib_disc/spatial_discretization/ip_data/ip_data.h"
#include "lib_disc/spatial_discretization/ip_data/user_function.h"
#include "lib_disc/spatial_discretization/ip_data/data_linker.h"

#include "lua_util.h"

using namespace std;

namespace ug
{
namespace bridge
{

///////////////////////////////////////////////////////////////////////////////
// Implementation
///////////////////////////////////////////////////////////////////////////////

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
	: public IPData<TData, dim>,
	  public boost::function<void (TData& res, const MathVector<dim>& x,number time)>
{
	///	Base class type
		typedef IPData<TData, dim> base_type;

	///	Functor type
		typedef boost::function<void (TData& res, const MathVector<dim>& x,number time)> func_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::value;

	public:
	///	Constructor
		LuaUserData() : func_type(boost::ref(*this))
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
				std::stringstream ss;
				ss << "ERROR in 'LuaUserData::operator(...)': Error while "
						"running callback '" << m_callbackName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n";
				throw(UGError(true, ss.str().c_str()));
			}

		//	read return value
			lua_traits<TData>::read(m_L, D);

		//	pop values
			lua_pop(m_L, retSize);
		}

	///	implement as a IPData
		virtual bool compute(bool bDeriv = false)
		{
			for(size_t s = 0; s < num_series(); ++s)
				for(size_t i = 0; i < num_ip(s); ++i)
				{
					this->operator()(	value(s,i),
										ip(s, i),
										time());
				}
			return true;
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
	: public boost::function<bool (TData& res, const MathVector<dim>& x,number time)>
{
	///	Functor type
		typedef boost::function<bool (TData& res, const MathVector<dim>& x,number time)> func_type;

	public:
	///	Constructor
		LuaBoundaryData()  : func_type(boost::ref(*this))
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
				std::stringstream ss;
				ss << "ERROR in 'LuaBoundaryData::operator(...)': Error while "
						"running callback '" << m_callbackName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n";
				throw(UGError(true, ss.str().c_str()));
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
	: public DataLinkerEqualData<TData, dim, TDataIn>
{
	public:
	//	type of base class
		typedef DataLinkerEqualData<TData, dim, TDataIn> base_type;

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
		using base_type::input_num_fct;
		using base_type::input_common_fct;

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
				std::stringstream ss;
				ss << "ERROR in 'LuaUserFunction::operator(...)': Error while "
						"running callback '" << m_cbValueName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n";
				throw(UGError(true, ss.str().c_str()));
			}

		//	read return value
			lua_traits<TData>::read(m_L, out);

		//	pop values
			lua_pop(m_L, retSize);
		}

	///	computes the value
		virtual bool compute(bool bDeriv)
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
			if(!bDeriv || this->zero_derivative()) return true;

		//	clear all derivative values
			this->clear_derivative_values();

		//	loop all inputs
			for(size_t c = 0; c < vDataIn.size(); ++c)
			{
			//	check if input has derivative
				if(this->zero_derivative(c)) continue;

			//	loop ips
				for(size_t s = 0; s < num_series(); ++s)
					for(size_t ip = 0; ip < num_ip(s); ++ip)
					{
					//	gather all input data for this ip
						vDataIn[c] = input_value(c, s, ip);

					//	data of derivative w.r.t. one component at ip-values
						TData derivVal;

					//	evaluate data at ip
						eval_deriv(derivVal, vDataIn, c);

					//	loop functions
						for(size_t fct = 0; fct < input_num_fct(c); ++fct)
						{
						//	get common fct id for this function
							const size_t commonFct = input_common_fct(c, fct);

						//	loop dofs
							for(size_t dof = 0; dof < num_sh(s, fct); ++dof)
							{
								linker_traits<TData, TDataIn>::
								mult_add(deriv(s, ip, commonFct, dof),
								         derivVal,
								         input_deriv(c, s, ip, fct, dof));
							}
						}
					}
			}

			return true;
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
				std::stringstream ss;
				ss << "ERROR in 'LuaUserFunction::eval_value(...)': Error while "
						"running callback '" << m_cbValueName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n";
				throw(UGError(true, ss.str().c_str()));
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
				std::stringstream ss;
				ss << "ERROR in 'LuaUserFunction::eval_deriv': Error while "
						"running callback '" << m_cbDerivName[arg] << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n";
				throw(UGError(true, ss.str().c_str()));
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
//		UG_LOG("Push value i=" << i << ": " << val<<"\n");
	}

	va_end(ap);


	if(lua_pcall(m_L, numArgs, 1, 0) != 0)
	{
		std::stringstream ss;
		ss << "ERROR in 'LuaUserNumberNumberFunction::operator(...)': Error while "
				"running callback '" << m_callbackName << "',"
				" lua message: "<< lua_tostring(m_L, -1) << "\n";
		throw(UGError(true, ss.str().c_str()));
	}

	number c = luaL_checknumber(m_L, -1);
	lua_pop(m_L, 1);

	return c;
}


////////////////////////////////////////////////////////////////////////////////
// LuaFunction
////////////////////////////////////////////////////////////////////////////////

/**
 * \tparam		TData		Return value type
 * \tparam		TDataIn		Input daten type
 */
template <typename TData, typename TDataIn>
class LuaFunction
	: public IFunction<TData, TDataIn>
{
	public:
	///	constructor
		LuaFunction() : m_numArgs(0)
		{
			m_L = ug::script::GetDefaultLuaState();
			m_cbValueRef = LUA_NOREF;
		}

	///	sets the Lua function used to compute the data
	/**
	 * This function sets the lua callback. The name of the function is
	 * passed as a string. Make sure, that the function name is defined
	 * when executing the script.
	 */
		void set_lua_callback(const char* luaCallback, size_t numArgs)
		{
		//	store name (string) of callback
			m_cbValueName = luaCallback;

		//	obtain a reference
			lua_getglobal(m_L, m_cbValueName);

		//	store reference to lua function
			m_cbValueRef = luaL_ref(m_L, LUA_REGISTRYINDEX);

		//	remember number of arguments to be used
			m_numArgs = numArgs;
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
				std::stringstream ss;
				ss << "ERROR in 'LuaUserFunction::operator(...)': Error while "
						"running callback '" << m_cbValueName << "',"
						" lua message: "<< lua_tostring(m_L, -1) << "\n";
				throw(UGError(true, ss.str().c_str()));
			}

		//	read return value
			lua_traits<TData>::read(m_L, out);

		//	pop values
			lua_pop(m_L, retSize);
		}

	protected:
	///	callback name as string
		const char* m_cbValueName;

	///	reference to lua function
		int m_cbValueRef;

	///	lua state
		lua_State*	m_L;

	///	number of arguments to use
		size_t m_numArgs;
};

///////////////////////////////////////////////////////////////////////////////
// Registration
///////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim>
bool RegisterLuaUserDataType(Registry& reg, string type, const char* parentGroup)
{
	string grp = std::string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

//	LuaUser"Type"
	{
		typedef LuaUserData<TData, dim> T;
		typedef IPData<TData, dim> TBase;
		typedef boost::function<void (TData& res, const MathVector<dim>& x,number time)> TBase2;
		string name = string("LuaUser").append(type).append(dimSuffix);
		reg.add_class_<T, TBase, TBase2>(name, grp)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
		reg.add_class_to_group(name, string("LuaUser").append(type), dimTag);
	}

//	LuaBoundary"Type"
	{
		typedef LuaBoundaryData<TData, dim> T;
		typedef boost::function<bool (TData& res, const MathVector<dim>& x,number time)> TBase;
		string name = string("LuaBoundary").append(type).append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
		reg.add_class_to_group(name, string("LuaBoundary").append(type), dimTag);
	}

	return true;
}

template <int dim>
void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{
	string grp = std::string(parentGroup);

	string dimSuffix = GetDomainSuffix<dim>();
	string dimTag = GetDomainTag<dim>();

	RegisterLuaUserDataType<number, dim>(reg, "Number", parentGroup);
	RegisterLuaUserDataType<MathVector<dim>, dim>(reg, "Vector", parentGroup);
	RegisterLuaUserDataType<MathMatrix<dim,dim>, dim>(reg, "Matrix", parentGroup);

//	LuaUserFunctionNumber
	{
		typedef LuaUserFunction<number, dim, number> T;
		typedef DataLinkerEqualData<number, dim, number> TBase;
		string name = string("LuaUserFunctionNumber").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_lua_value_callback", &T::set_lua_value_callback)
			.add_method("set_lua_deriv_callback", &T::set_lua_deriv_callback);
		reg.add_class_to_group(name, "LuaUserFunctionNumber", dimTag);
	}

//	LuaUserFunctionMatrixNumber
	{
		typedef LuaUserFunction<MathMatrix<dim,dim>, dim, number> T;
		typedef DataLinkerEqualData<MathMatrix<dim,dim>, dim, number> TBase;
		string name = string("LuaUserFunctionMatrixNumber").append(dimSuffix);
		reg.add_class_<T, TBase>(name, grp)
			.add_constructor()
			.add_method("set_lua_value_callback", &T::set_lua_value_callback)
			.add_method("set_lua_deriv_callback", &T::set_lua_deriv_callback);
		reg.add_class_to_group(name, "LuaUserFunctionMatrixNumber", dimTag);
	}

}

void RegisterLuaUserData(Registry& reg, const char* parentGroup)
{

//	LuaUserNumberNumberFunction
	{
		typedef LuaUserNumberNumberFunction T;
		reg.add_class_<T>("LuaUserNumberNumberFunction", parentGroup)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

//	LuaFunctionNumber
	{
		typedef LuaFunction<number, number> T;
		typedef IFunction<number, number> TBase;
		reg.add_class_<T, TBase>("LuaFunctionNumber", parentGroup)
			.add_constructor()
			.add_method("set_lua_callback", &T::set_lua_callback);
	}

#ifdef UG_DIM_1
	RegisterLuaUserData<1>(reg, parentGroup);
#endif
#ifdef UG_DIM_2
	RegisterLuaUserData<2>(reg, parentGroup);
#endif
#ifdef UG_DIM_3
	RegisterLuaUserData<3>(reg, parentGroup);
#endif
}

} // end namespace
} // end namepace

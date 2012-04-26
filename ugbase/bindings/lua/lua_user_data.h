/*
 * lua_user_data.h
 *      Author: andreasvogel
 */

#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__

#include <stdarg.h>
#include "registry/registry.h"

extern "C" {
#include "externals/lua/lua.h"
}

#include "lua_util.h"

#include "common/common.h"
#include "common/math/ugmath.h"
#include "lib_disc/spatial_disc/ip_data/ip_data.h"
#include "lib_disc/spatial_disc/ip_data/user_function.h"
#include "lib_disc/spatial_disc/ip_data/data_linker.h"

namespace ug
{

////////////////////////////////////////////////////////////////////////////////
// LuaUserData
////////////////////////////////////////////////////////////////////////////////

///	returns true if callback exists
bool CheckLuaCallbackName(const char* name);

// predeclaration
template <typename TData, int dim, typename TRet = void>
class LuaUserDataFactory;

/// provides data specified in the lua script
/**
 * This class implements the IPData interface to provide data at arbitrary
 * points. A Lua callback is used to evaluate the data.
 *
 * NOTE: If the LuaUserData has been created by the LuaUserDataFactory, then
 * 		 the fromFactory flag is set to true internally and while deconstruction
 * 		 of the instance LuaUserDataFactory::remove is called.
 */
template <typename TData, int dim, typename TRet = void>
class LuaUserData : public IPData<TData, dim, TRet>
{
	///	friend class
		friend class LuaUserDataFactory<TData, dim, TRet>;

	///	Base class type
		typedef IPData<TData, dim, TRet> base_type;

		using base_type::num_series;
		using base_type::num_ip;
		using base_type::ip;
		using base_type::time;
		using base_type::subset;
		using base_type::value;

	public:
	///	Constructor
	/**
	 * Creates a LuaUserData that uses a Lua function to evaluate some data.
	 * NOTE: The Lua callback function is called once with dummy parameters
	 * 		 in order to check the correct return values.
	 *
	 * @param luaCallback		Name of Lua Callback Function
	 */
		LuaUserData(const char* luaCallback);

	///	destructor: frees lua callback, unregisters from LuaUserDataFactory if used
		virtual ~LuaUserData();

	///	returns string of required callback signature
		static std::string signature();

	///	returns name of UserData
		static std::string name();

	///	returns true if callback has correct return values
		static bool check_callback_returns(const char* callName,
										   const bool bThrow = false);

	///	evaluates the data at a given point and time
		TRet operator() (TData& D, const MathVector<dim>& x, number time, int si) const;

	///	implement as a IPData
		virtual void compute(bool bDeriv = false);

	protected:
	///	sets that LuaUserData is created by LuaUserDataFactory
		void set_created_from_factory(bool bFromFactory) {m_bFromFactory = bFromFactory;}

	protected:
	///	callback name as string
		std::string m_callbackName;

	///	reference to lua function
		int m_callbackRef;

	///	flag, indicating if created from factory
		bool m_bFromFactory;

	///	lua state
		lua_State*	m_L;
};

////////////////////////////////////////////////////////////////////////////////
// LuaUserDataFactory
////////////////////////////////////////////////////////////////////////////////

/// Factory providing LuaUserData
/**
 * This class is a singleton providing LuaUserData for a given callback name.
 * The factory creates a new LuaUserData and returns it as a SmartPtr if no
 * Data has been created with the same callback name before. It already a
 * LuaUserData using the same callback exists, then a SmartPtr to that instance
 * is returned. Internally, therefore plain pointers and reference counters
 * of the created SmartPtrs are stored and LuaUserData unregister from the
 * List of created LuaUserData in this factory when their destructor is called.
 *
 * \tparam	TData	Data produced by LuaUserData
 * \tparam	dim		world dimension
 * \tparam	TRet	bool flag or void
 */
template <typename TData, int dim, typename TRet>
class LuaUserDataFactory
{
	friend class LuaUserData<TData,dim,TRet>;

	protected:
	///	private constructor, since singleton
		LuaUserDataFactory(){};

	//	disallow copy
		LuaUserDataFactory(const LuaUserDataFactory&);
		LuaUserDataFactory& operator=(const LuaUserDataFactory&);

	///	singleton provider
		static LuaUserDataFactory<TData,dim,TRet>& instance()
		{
			static LuaUserDataFactory<TData,dim,TRet> inst;
			return inst;
		}

	///	returns new Data if not already created, already existing else
		static SmartPtr<LuaUserData<TData,dim,TRet> > provide_or_create(std::string name);

	///	removes the user data
		static void remove(std::string name);

	///	storage of already created data
		static std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> > m_mData;

	public:
	/**
	 * Returns a LuaUserDataas a SmartPtr for a callback name. If a LuaUserData has
	 * already been created by the factory it simply returns a new SmartPtr
	 * without creation, else the LuaUserData is created. All LuaUserData
	 * store internally, if they have been created by a factory, and if so
	 * they unregister from the factory when their destructor is called.
	 *
	 * NOTE: In order to allow any garbage collection of the LuaUserData,
	 * therefore this factory can not store the LuaUserData as a SmartPtr, but
	 * as a plain pointer together with the reference counter.
	 *
	 * @param name		lua callback name
	 * @return			LuaUserNumber as SmartPtr
	 */
		static SmartPtr<LuaUserData<TData,dim,TRet> > create(std::string name)
		{
			instance().provide_or_create(name);
		}
};

////////////////////////////////////////////////////////////////////////////////
// LuaUserFunction
////////////////////////////////////////////////////////////////////////////////

/// maps several data values to an output data value using a lua callback
/**
 * This class provides the evaluation of a user function, that is specified
 * in the script. Several data (of the same c++-type) can be used as input,
 * a data (of the same type) is returned.
 */
template <typename TData, int dim, typename TDataIn>
class LuaUserFunction : public DataLinkerEqualData<TData, dim, TDataIn>
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
		LuaUserFunction(const char* luaCallback, size_t numArgs);

	///	destructor frees the reference
		virtual ~LuaUserFunction();

	///	sets the Lua function used to compute the derivative
		void set_deriv(size_t arg, const char* luaCallback);

	///	evaluates the data
		virtual void operator() (TData& out, int numArgs, ...);

	///	computes the value
		virtual void compute(bool bDeriv);

	protected:
	///	sets the Lua function used to compute the data
		void set_lua_value_callback(const char* luaCallback, size_t numArgs);

	///	frees the callback-reference, if a callback was set.
		void free_callback_ref();

	///	frees callback-references for derivate callbacks
		void free_deriv_callback_ref(size_t arg);

	///	evaluates the data
		void eval_value(TData& out, const std::vector<TDataIn>& dataIn);

	///	evaluates the data
		void eval_deriv(TData& out, const std::vector<TDataIn>& dataIn, size_t arg);

	protected:
	///	callback name as string
		std::string m_cbValueName;
		std::vector<std::string> m_cbDerivName;

	///	reference to lua function
		int m_cbValueRef;
		std::vector<int> m_cbDerivRef;

	///	lua state
		lua_State*	m_L;

	///	number of arguments to use
		size_t m_numArgs;
};


/// this class maps a scalar value an output scalar value using a lua callback
class LuaUserNumberNumberFunction
{
	public:
		LuaUserNumberNumberFunction();

		void set_lua_callback(const char* luaCallback);

		number operator() ( int numArgs, ... ) const;

	protected:
		std::string m_callbackName;
		int m_callbackRef;
		lua_State*	m_L;
};


////////////////////////////////////////////////////////////////////////////////
// LuaFunction
////////////////////////////////////////////////////////////////////////////////

/**
 * \tparam		TData		Return value type
 * \tparam		TDataIn		Input daten type
 */
template <typename TData, typename TDataIn>
class LuaFunction : public IFunction<TData, TDataIn>
{
	public:
	///	constructor
		LuaFunction();

	///	sets the Lua function used to compute the data
	/**
	 * This function sets the lua callback. The name of the function is
	 * passed as a string. Make sure, that the function name is defined
	 * when executing the script.
	 */
		void set_lua_callback(const char* luaCallback, size_t numArgs);

	///	evaluates the data
		virtual void operator() (TData& out, int numArgs, ...);

	protected:
	///	callback name as string
		std::string m_cbValueName;

	///	reference to lua function
		int m_cbValueRef;

	///	lua state
		lua_State*	m_L;

	///	number of arguments to use
		size_t m_numArgs;
};



namespace bridge{

void RegisterLuaUserData(Registry& reg, const char* parentGroup);
void RegisterLuaCondUserNumber(Registry& reg, const char* parentGroup);

} // end namepace bridge
} // end namespace ug

#include "lua_user_data_impl.h"

#endif /* __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA__ */

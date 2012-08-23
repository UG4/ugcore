/*
 * lua_user_data_impl.h
 *
 *  Created on: 26.04.2012
 *      Author: andreasvogel
 */

#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA_IMPL_
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA_IMPL_

#include "lua_user_data.h"
#include "lib_disc/spatial_disc/user_data/data_linker_traits.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"

namespace ug{

///////////////////////////////////////////////////////////////////////////////
// Lua Traits
///////////////////////////////////////////////////////////////////////////////

/// Helper to access a return value on the stack.
/**	If the value can't be converted to a number, an error is thrown*/
inline number ReturnValueToNumber(lua_State* L, int index){
	if(!lua_isnumber(L, index))
		UG_THROW("ReturnValueToNumber: Data passed from Lua: "
						"Can't convert return value to number!");
	return lua_tonumber(L, index);
}

/// Helper to access a return value on the stack.
/**	If the value can't be converted to a boolean, an error is thrown*/
inline number ReturnValueToBool(lua_State* L, int index){
	if(!lua_isboolean(L, index))
		UG_THROW("ReturnValueToBool: Data passed from Lua: "
						"Can't convert return value to boolean!");

	return lua_toboolean(L, index);
}

/// Helper to access a return value on the stack.
/**	If the value can't be converted to a integer, an error is thrown*/
inline int ReturnValueToInteger(lua_State* L, int index){
	if(!lua_isnumber(L, index))
		UG_THROW("ReturnValueToBool: Data passed from Lua: "
						"Can't convert return value to integer!");

	return lua_tointeger(L, index);
}



/// Lua Traits to push/pop on lua stack
template <typename TData>
struct lua_traits;

template <>
struct lua_traits<void>
{
	static const int size = 0;

	static void push(lua_State*	L, const bool&){}

	static void read(lua_State* L, bool&, int index = -1){}

	static void do_return(const bool&) {return;}

	static bool check(lua_State* L, int index = -1){return true;}

	static std::string signature()
	{
		return "void";
	}

	static std::string name()
	{
		return "void";
	}
};

template <>
struct lua_traits<bool>
{
	static const int size = 1;

	static void push(lua_State*	L, const bool& b)
	{
		lua_pushboolean(L, b);
	}

	static void read(lua_State* L, bool& b, int index = -1)
	{
		b =  ReturnValueToBool(L, index);
	}

	static bool check(lua_State* L, int index = -1)
	{
		return lua_isboolean(L, index);
	}

	static bool do_return(const bool& b) {return b;}

	static std::string signature()
	{
		return "bool";
	}

	static std::string name()
	{
		return "Bool";
	}
};

template <>
struct lua_traits<int>
{
	static const int size = 1;

	static void push(lua_State*	L, const int& c)
	{
		lua_pushinteger(L, c);
	}

	static void read(lua_State* L, int& c, int index = -1)
	{
		c = ReturnValueToInteger(L, index);
	}

	static bool check(lua_State* L, int index = -1)
	{
		return lua_isnumber(L, index);
	}

	static std::string signature()
	{
		return "integer";
	}

	static std::string name()
	{
		return "Integer";
	}
};

template <>
struct lua_traits<number>
{
	static const int size = 1;

	static void push(lua_State*	L, const number& c)
	{
		lua_pushnumber(L, c);
	}

	static void read(lua_State* L, number& c, int index = -1)
	{
		c = ReturnValueToNumber(L, index);
	}

	static bool check(lua_State* L, int index = -1)
	{
		return lua_isnumber(L, index);
	}

	static std::string signature()
	{
		return "number";
	}

	static std::string name()
	{
		return "Number";
	}
};

template <std::size_t dim>
struct lua_traits< ug::MathVector<dim> >
{
	static const int size = dim;

	static void push(lua_State*	L, const MathVector<dim>& x)
	{
		for(size_t i = 0; i < dim; ++i)
			lua_pushnumber(L, x[i]);
	}

	static void read(lua_State* L, MathVector<dim>& x, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
				x[dim-1-i] = ReturnValueToNumber(L, index--);
		}
	}

	static bool check(lua_State* L, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
			if(!lua_isnumber(L, index--)) return false;
		}
		return true;
	}

	static std::string signature()
	{
		static char cmp[] = {'x', 'y', 'z'};
		std::stringstream ss;
		for(size_t i = 0; i < dim; ++i){
			if(i != 0) ss << ", ";
			ss << "v" << cmp[i];
		}
		return ss.str();
	}

	static std::string name()
	{
		std::stringstream ss;
		ss << "Vector";
		return ss.str();
	}
};


template <std::size_t dim>
struct lua_traits< MathMatrix<dim, dim> >
{
	static const int size = dim*dim;

	static void push(lua_State*	L, const MathMatrix<dim, dim>& D)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				lua_pushnumber(L, D[i][j]);
			}
		}

	}

	static void read(lua_State* L, MathMatrix<dim, dim>& D, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				D[dim-1-j][dim-1-i] = ReturnValueToNumber(L, index--);
			}
		}
	}

	static bool check(lua_State* L, int index = -1)
	{
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				if(!lua_isnumber(L, index--)) return false;
			}
		}
		return true;
	}

	static std::string signature()
	{
		static char cmp[] = {'x', 'y', 'z'};
		std::stringstream ss;
		for(size_t i = 0; i < dim; ++i){
			for(size_t j = 0; j < dim; ++j){
				if(i != 0 || j != 0) ss << ", ";
				ss << "D" << cmp[i] << cmp[j];
			}
		}
		return ss.str();
	}

	static std::string name()
	{
		std::stringstream ss;
		ss << "Matrix";
		return ss.str();
	}
};


////////////////////////////////////////////////////////////////////////////////
// LuaUserData
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TRet>
std::string LuaUserData<TData,dim,TRet>::signature()
{
	std::stringstream ss;
	ss << "function name(";
	if(dim >= 1) ss << "x";
	if(dim >= 2) ss << ", y";
	if(dim >= 3) ss << ", z";
	ss << ", t, si)\n   ... \n   return ";
	if(lua_traits<TRet>::size != 0)
		ss << lua_traits<TRet>::signature() << ", ";
	ss << lua_traits<TData>::signature();
	ss << "\nend";
	return ss.str();
}


template <typename TData, int dim, typename TRet>
std::string LuaUserData<TData,dim,TRet>::name()
{
	std::stringstream ss;
	ss << "Lua";
	if(lua_traits<TRet>::size > 0) ss << "Cond";
	ss << "User" << lua_traits<TData>::name() << dim << "d";
	return ss.str();
}

template <typename TData, int dim, typename TRet>
LuaUserData<TData,dim,TRet>::LuaUserData(const char* luaCallback)
	: m_callbackName(luaCallback), m_bFromFactory(false)
{
//	make a test run
	check_callback_returns(luaCallback, true);

//	get lua state
	m_L = ug::script::GetDefaultLuaState();

//	obtain a reference
	lua_getglobal(m_L, m_callbackName.c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW(name() << ": Specified lua callback "
						"does not exist: " << m_callbackName);
	}

//	store reference to lua function
	m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);
}

template <typename TData, int dim, typename TRet>
bool LuaUserData<TData,dim,TRet>::
check_callback_returns(const char* callName, const bool bThrow)
{
//	get lua state
	lua_State* L = ug::script::GetDefaultLuaState();

//	get current stack level
	const int level = lua_gettop(L);

//	obtain a reference
	lua_getglobal(L, callName);

//	check if reference is valid
	if(lua_isnil(L, -1)) {
		if(bThrow) {
			UG_THROW(name() << ": Cannot find specified lua callback "
							" with name: "<<callName);
		}
		else {
			return false;
		}
	}

//	get reference
	int callbackRef = luaL_ref(L, LUA_REGISTRYINDEX);

//	dummy values to invoke the callback once
	MathVector<dim> x; x = 0.0;
	number time = 0.0;
	int si = 0;

//	push the callback function on the stack
	lua_rawgeti(L, LUA_REGISTRYINDEX, callbackRef);

//  push space coordinates on stack
	lua_traits<MathVector<dim> >::push(L, x);

//	push time on stack
	lua_traits<number>::push(L, time);

//	push subset on stack
	lua_traits<int>::push(L, si);

//	compute total args size
	const int argSize = lua_traits<MathVector<dim> >::size
							+ lua_traits<number>::size
							+ lua_traits<int>::size;

//	compute total return size
	const int retSize = lua_traits<TData>::size + lua_traits<TRet>::size;

//	call lua function
	if(lua_pcall(L, argSize, LUA_MULTRET, 0) != 0)
		UG_THROW(name() << ": Error while "
						"testing callback '" << callName << "',"
						" lua message: "<< lua_tostring(L, -1));

//	get number of results
	const int numResults = lua_gettop(L) - level;

//	success flag
	bool bRet = true;

//	if number of results is wrong return error
	if(numResults != retSize){
		if(bThrow){
			UG_THROW(name() << ": Number of return values incorrect "
							"for callback '"<<callName<<"'. "
							"Required: "<<retSize<<", passed: "<<numResults
							<<". Use signature as follows:\n"
							<< signature());
		}
		else{
			bRet = false;
		}
	}

//	check return value
	if(!lua_traits<TData>::check(L)){
		if(bThrow){
			UG_THROW(name() << ": Data values type incorrect "
							"for callback '"<<callName<<"'. "
							"Use signature as follows:\n"
							<< signature());
		}
		else{
			bRet = false;
		}
	}

//	read return flag (may be void)
	if(!lua_traits<TRet>::check(L, -retSize)){
		if(bThrow){
			UG_THROW("LuaUserData: Return values type incorrect "
							"for callback '"<<callName<<"'. "
							"Use signature as follows:\n"
							<< signature());
		}
		else{
			bRet = false;
		}
	}

//	pop values
	lua_pop(L, numResults);

//	free reference to callback
	luaL_unref(L, LUA_REGISTRYINDEX, callbackRef);

//	return match
	return bRet;
}

template <typename TData, int dim, typename TRet>
TRet LuaUserData<TData,dim,TRet>::
evaluate(TData& D, const MathVector<dim>& x, number time, int si) const
{
//	push the callback function on the stack
	lua_rawgeti(m_L, LUA_REGISTRYINDEX, m_callbackRef);

//  push space coordinates on stack
	lua_traits<MathVector<dim> >::push(m_L, x);

//	push time on stack
	lua_traits<number>::push(m_L, time);

//	push subset index on stack
	lua_traits<int>::push(m_L, si);

//	compute total args size
	const int argSize = lua_traits<MathVector<dim> >::size
						+ lua_traits<number>::size
						+ lua_traits<int>::size;

//	compute total return size
	const int retSize = lua_traits<TData>::size + lua_traits<TRet>::size;

//	call lua function
	if(lua_pcall(m_L, argSize, retSize, 0) != 0)
		UG_THROW(name() << "::operator(...): Error while "
						"running callback '" << m_callbackName << "',"
						" lua message: "<< lua_tostring(m_L, -1)<<".\n"
						"Use signature as follows:\n"
						<< signature());

	bool res = false;
	try{
	//	read return value
		lua_traits<TData>::read(m_L, D);

	//	read return flag (may be void)
		lua_traits<TRet>::read(m_L, res, -retSize);
	}
	UG_CATCH_THROW(name() << "::operator(...): Error while running "
					"callback '" << m_callbackName << "'."
					"Use signature as follows:\n"
					<< signature());

//	pop values
	lua_pop(m_L, retSize);

//	forward flag
	return lua_traits<TRet>::do_return(res);
}

template <typename TData, int dim, typename TRet>
LuaUserData<TData,dim,TRet>::~LuaUserData()
{
//	free reference to callback
	luaL_unref(m_L, LUA_REGISTRYINDEX, m_callbackRef);

	if(m_bFromFactory)
		LuaUserDataFactory<TData,dim,TRet>::remove(m_callbackName);
}

////////////////////////////////////////////////////////////////////////////////
// LuaUserDataFactory
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TRet>
SmartPtr<LuaUserData<TData,dim,TRet> >
LuaUserDataFactory<TData,dim,TRet>::provide_or_create(std::string name)
{
	typedef std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> > Map;
	typedef typename Map::iterator iterator;

//	check for element
	iterator iter = m_mData.find(name);

//	if name does not exist, create new one
	if(iter == m_mData.end())
	{
		SmartPtr<LuaUserData<TData,dim,TRet> > sp
			= CreateSmartPtr(new LuaUserData<TData,dim,TRet>(name.c_str()));

	//	the LuaUserData must remember to unregister itself at destruction
		sp->set_created_from_factory(true);

	//	NOTE AND WARNING: This is very hacky and dangerous. We only do this
	//	since we exactly know what we are doing  and everything is save and
	//	only in protected or private area. However, if you once want to change
	//	this code, please be aware, that we store here plain pointers and
	//	associated reference counters of a SmartPtr. This should not be done
	//	in general and this kind of coding is not recommended at all. Please
	//	use different approaches whenever possible.
		std::pair<LuaUserData<TData,dim,TRet>*, int*>& data = m_mData[name];
		data.first = sp.get();
		data.second = sp.refcount_ptr();

		return sp;
	}
//	else return present data
	{
	//	NOTE AND WARNING: This is very hacky and dangerous. We only do this
	//	since we exactly know what we are doing  and everything is save and
	//	only in protected or private area. However, if you once want to change
	//	this code, please be aware, that we store here plain pointers and
	//	associated reference counters of a SmartPtr. This should not be done
	//	in general and this kind of coding is not recommended at all. Please
	//	use different approaches whenever possible.
		std::pair<LuaUserData<TData,dim,TRet>*, int*>& data = iter->second;
		return SmartPtr<LuaUserData<TData,dim,TRet> >(data.first, data.second);
	}
}

template <typename TData, int dim, typename TRet>
void
LuaUserDataFactory<TData,dim,TRet>::remove(std::string name)
{
	typedef std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> > Map;
	typedef typename Map::iterator iterator;

//	check for element
	iterator iter = m_mData.find(name);

//	if name does not exist, create new one
	if(iter == m_mData.end())
		UG_THROW("LuaUserDataFactory: trying to remove non-registered"
						" data with name: "<<name);

	m_mData.erase(iter);
}


// instantiation of static member
template <typename TData, int dim, typename TRet>
std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> >
LuaUserDataFactory<TData,dim,TRet>::m_mData = std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> >();

////////////////////////////////////////////////////////////////////////////////
// LuaUserFunction
////////////////////////////////////////////////////////////////////////////////

template <typename TData, int dim, typename TDataIn>
LuaUserFunction<TData,dim,TDataIn>::
LuaUserFunction(const char* luaCallback, size_t numArgs)
	: m_numArgs(numArgs), m_bPosTimeNeed(false)
{
	m_L = ug::script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
	m_cbDerivRef.clear();
	m_cbDerivName.clear();
	set_lua_value_callback(luaCallback, numArgs);
}

template <typename TData, int dim, typename TDataIn>
LuaUserFunction<TData,dim,TDataIn>::
LuaUserFunction(const char* luaCallback, size_t numArgs, bool bPosTimeNeed)
	: m_numArgs(numArgs), m_bPosTimeNeed(bPosTimeNeed)
{
	m_L = ug::script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
	m_cbDerivRef.clear();
	m_cbDerivName.clear();
	set_lua_value_callback(luaCallback, numArgs);
}

template <typename TData, int dim, typename TDataIn>
LuaUserFunction<TData,dim,TDataIn>::~LuaUserFunction()
{
//	free reference to callback
	free_callback_ref();

//	free references to derivate callbacks
	for(size_t i = 0; i < m_numArgs; ++i){
		free_deriv_callback_ref(i);
	}
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::free_callback_ref()
{
	if(m_cbValueRef != LUA_NOREF){
		luaL_unref(m_L, LUA_REGISTRYINDEX, m_cbValueRef);
		m_cbValueRef = LUA_NOREF;
	}
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::free_deriv_callback_ref(size_t arg)
{
	if(m_cbDerivRef[arg] != LUA_NOREF){
		luaL_unref(m_L, LUA_REGISTRYINDEX, m_cbDerivRef[arg]);
		m_cbDerivRef[arg] = LUA_NOREF;
	}
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::set_lua_value_callback(const char* luaCallback, size_t numArgs)
{
//	store name (string) of callback
	m_cbValueName = luaCallback;

//	obtain a reference
	lua_getglobal(m_L, m_cbValueName.c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW("LuaUserFunction::set_lua_value_callback(...):"
				"Specified callback does not exist: " << m_cbValueName);
	}

//	if a callback was already set, we have to free the old one
	free_callback_ref();

//	store reference to lua function
	m_cbValueRef = luaL_ref(m_L, LUA_REGISTRYINDEX);

//	remember number of arguments to be used
	m_numArgs = numArgs;
	m_cbDerivName.resize(numArgs);
	m_cbDerivRef.resize(numArgs, LUA_NOREF);

//	set num inputs for linker
	set_num_input(numArgs);
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::set_deriv(size_t arg, const char* luaCallback)
{
//	check number of arg
	if(arg >= m_numArgs)
		UG_THROW("LuaUserFunction::set_lua_deriv_callback: Trying "
				"to set a derivative for argument " << arg <<", that "
				"does not exist. Number of arguments is "<<m_numArgs);

//	store name (string) of callback
	m_cbDerivName[arg] = luaCallback;

//	free old reference
	free_deriv_callback_ref(arg);

//	obtain a reference
	lua_getglobal(m_L, m_cbDerivName[arg].c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW("LuaUserFunction::set_lua_deriv_callback(...):"
				"Specified callback does not exist: " << m_cbDerivName[arg]);
	}

//	store reference to lua function
	m_cbDerivRef[arg] = luaL_ref(m_L, LUA_REGISTRYINDEX);
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::operator() (TData& out, int numArgs, ...) const
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
		UG_THROW("LuaUserFunction::operator(...): Error while "
					"running callback '" << m_cbValueName << "',"
					" lua message: "<< lua_tostring(m_L, -1));

	try{
	//	read return value
		lua_traits<TData>::read(m_L, out);
	}
	UG_CATCH_THROW("LuaUserFunction::operator(...): Error while running "
					"callback '" << m_cbValueName << "'");

//	pop values
	lua_pop(m_L, retSize);
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::eval_value(TData& out, const std::vector<TDataIn>& dataIn,
													const MathVector<dim>& x, number time, int si) const
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

//	if needed, read additional coordinate, time and subset index arguments and push them to the lua stack
	if(m_bPosTimeNeed){
		lua_traits<MathVector<dim> >::push(m_L, x);
		lua_traits<number>::push(m_L, time);
		lua_traits<int>::push(m_L, si);
	}

//	compute total args size
	size_t argSize = lua_traits<TDataIn>::size * dataIn.size();
	if(m_bPosTimeNeed){
		argSize += 	lua_traits<MathVector<dim> >::size
					+ lua_traits<number>::size
					+ lua_traits<int>::size;
	}

//	compute total return size
	size_t retSize = lua_traits<TData>::size;

//	call lua function
	if(lua_pcall(m_L, argSize, retSize, 0) != 0)
		UG_THROW("LuaUserFunction::eval_value(...): Error while "
					"running callback '" << m_cbValueName << "',"
					" lua message: "<< lua_tostring(m_L, -1));

	try{
	//	read return value
		lua_traits<TData>::read(m_L, out);
	}
	UG_CATCH_THROW("LuaUserFunction::eval_value(...): Error while "
					"running callback '" << m_cbValueName << "'");

//	pop values
	lua_pop(m_L, retSize);
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::eval_deriv(TData& out, const std::vector<TDataIn>& dataIn,
		 	 	 	 	 	 	 	 	 	 	 	 const MathVector<dim>& x, number time, int si, size_t arg) const
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

//	if needed, read additional coordinate, time and subset index arguments and push them to the lua stack
	if(m_bPosTimeNeed){
		lua_traits<MathVector<dim> >::push(m_L, x);
		lua_traits<number>::push(m_L, time);
		lua_traits<int>::push(m_L, si);
	}

//	compute total args size
	size_t argSize = lua_traits<TDataIn>::size * dataIn.size();
	if(m_bPosTimeNeed){
		argSize += 	lua_traits<MathVector<dim> >::size
					+ lua_traits<number>::size
					+ lua_traits<int>::size;
	}

//	compute total return size
	size_t retSize = lua_traits<TData>::size;

//	call lua function
	if(lua_pcall(m_L, argSize, retSize, 0) != 0)
		UG_THROW("LuaUserFunction::eval_deriv: Error while "
				"running callback '" << m_cbDerivName[arg] << "',"
				" lua message: "<< lua_tostring(m_L, -1) );

	try{
	//	read return value
		lua_traits<TData>::read(m_L, out);
	}
	UG_CATCH_THROW("LuaUserFunction::eval_deriv(...): Error while "
				"running callback '" << m_cbDerivName[arg] << "'");

//	pop values
	lua_pop(m_L, retSize);
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::
evaluate (TData& value,
          const MathVector<dim>& globIP,
          number time, int si) const
{
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t c = 0; c < vDataIn.size(); ++c)
		(*m_vpUserData[c])(vDataIn[c], globIP, time, si);

//	evaluate data at ip
	eval_value(value, vDataIn, globIP, time, si);
}

template <typename TData, int dim, typename TDataIn>
template <int refDim>
void LuaUserFunction<TData,dim,TDataIn>::
evaluate (TData& value,
          const MathVector<dim>& globIP,
          number time, int si,
          LocalVector& u,
          GeometricObject* elem,
          const MathVector<dim> vCornerCoords[],
          const MathVector<refDim>& locIP) const
{
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t c = 0; c < vDataIn.size(); ++c)
		(*m_vpUserData[c])(vDataIn[c], globIP, time, si, u, elem, vCornerCoords, locIP);

//	evaluate data at ip
	eval_value(value, vDataIn, globIP, time, si);
}


template <typename TData, int dim, typename TDataIn>
template <int refDim>
void LuaUserFunction<TData,dim,TDataIn>::
evaluate(TData vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         LocalVector& u,
         GeometricObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<refDim> vLocIP[],
         const size_t nip,
         const MathMatrix<refDim, dim>* vJT) const
{
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t ip = 0; ip < nip; ++ip)
	{
		for(size_t c = 0; c < vDataIn.size(); ++c)
			(*m_vpUserData[c])(vDataIn[c], vGlobIP[ip], time, si, u, elem, vCornerCoords, vLocIP[ip]);

	//	evaluate data at ip
		eval_value(vValue[ip], vDataIn, vGlobIP[ip], time, si);
	}
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::
compute(LocalVector* u, GeometricObject* elem, bool bDeriv)
{
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

	const number t = this->time();
	const int si = this->subset();

	for(size_t s = 0; s < this->num_series(); ++s)
		for(size_t ip = 0; ip < this->num_ip(s); ++ip)
		{
		//	gather all input data for this ip
			for(size_t c = 0; c < vDataIn.size(); ++c)
				vDataIn[c] = m_vpUserData[c]->value(this->series_id(c,s), ip);

		//	evaluate data at ip
			eval_value(this->value(s,ip), vDataIn, this->ip(s, ip), t, si);
		}

//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return;

//	clear all derivative values
	this->clear_derivative_values();

//	loop all inputs
	for(size_t c = 0; c < vDataIn.size(); ++c)
	{
	//	check if input has derivative
		if(this->zero_derivative(c)) continue;

	//	loop ips
		for(size_t s = 0; s < this->num_series(); ++s)
			for(size_t ip = 0; ip < this->num_ip(s); ++ip)
			{
			//	gather all input data for this ip
				vDataIn[c] = m_vpUserData[c]->value(this->series_id(c,s), ip);

			//	data of derivative w.r.t. one component at ip-values
				TData derivVal;

			//	evaluate data at ip
				eval_deriv(derivVal, vDataIn, this->ip(s, ip), t, si, c);

			//	loop functions
				for(size_t fct = 0; fct < this->input_num_fct(c); ++fct)
				{
				//	get common fct id for this function
					const size_t commonFct = this->input_common_fct(c, fct);

				//	loop dofs
					for(size_t dof = 0; dof < this->num_sh(fct); ++dof)
					{
						linker_traits<TData, TDataIn>::
						mult_add(this->deriv(s, ip, commonFct, dof),
						         derivVal,
						         m_vpDependData[c]->deriv(this->series_id(c,s), ip, fct, dof));
					}
				}
			}
	}
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::set_num_input(size_t num)
{
//	resize arrays
	m_vpUserData.resize(num, NULL);
	m_vpDependData.resize(num, NULL);

//	forward size to base class
	base_type::set_num_input(num);
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::
set_input(size_t i, SmartPtr<UserData<TDataIn, dim> > data)
{
	UG_ASSERT(i < m_vpUserData.size(), "Input not needed");
	UG_ASSERT(i < m_vpDependData.size(), "Input not needed");

//	check input number
	if(i >= this->num_input())
		UG_THROW("DataLinker::set_input: Only " << this->num_input()
						<< " inputs can be set. Use 'set_num_input' to increase"
						" the number of needed inputs.");

//	remember userdata
	m_vpUserData[i] = data;

//	cast to dependent data
	m_vpDependData[i] = data.template cast_dynamic<DependentUserData<TDataIn, dim> >();

//	forward to base class
	base_type::set_input(i, data);
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::set_input(size_t i, number val)
{
	set_input(i, CreateConstUserData<dim>(val, TDataIn()));
}


////////////////////////////////////////////////////////////////////////////////
// LuaFunction
////////////////////////////////////////////////////////////////////////////////

template <typename TData, typename TDataIn>
LuaFunction<TData,TDataIn>::LuaFunction() : m_numArgs(0)
{
	m_L = ug::script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
}

template <typename TData, typename TDataIn>
void LuaFunction<TData,TDataIn>::set_lua_callback(const char* luaCallback, size_t numArgs)
{
//	store name (string) of callback
	m_cbValueName = luaCallback;

//	obtain a reference
	lua_getglobal(m_L, m_cbValueName.c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW("LuaFunction::set_lua_callback(...):"
				"Specified lua callback does not exist: " << m_cbValueName);
	}

//	store reference to lua function
	m_cbValueRef = luaL_ref(m_L, LUA_REGISTRYINDEX);

//	remember number of arguments to be used
	m_numArgs = numArgs;
}

template <typename TData, typename TDataIn>
void LuaFunction<TData,TDataIn>::operator() (TData& out, int numArgs, ...)
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
		UG_THROW("LuaFunction::operator(...): Error while "
					"running callback '" << m_cbValueName << "',"
					" lua message: "<< lua_tostring(m_L, -1));

	try{
	//	read return value
		lua_traits<TData>::read(m_L, out);
	}
	UG_CATCH_THROW("LuaFunction::operator(...): Error while running "
					"callback '" << m_cbValueName << "'");

//	pop values
	lua_pop(m_L, retSize);
}



} // end namespace ug

#endif /* LUA_USER_DATA_IMPL_H_ */

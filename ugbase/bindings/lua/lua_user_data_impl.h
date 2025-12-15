/*
 * Copyright (c) 2012-2015:  G-CSC, Goethe University Frankfurt
 * Author: Andreas Vogel
 * 
 * This file is part of UG4.
 * 
 * UG4 is free software: you can redistribute it and/or modify it under the
 * terms of the GNU Lesser General Public License version 3 (as published by the
 * Free Software Foundation) with the following additional attribution
 * requirements (according to LGPL/GPL v3 §7):
 * 
 * (1) The following notice must be displayed in the Appropriate Legal Notices
 * of covered and combined works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (2) The following notice must be displayed at a prominent place in the
 * terminal output of covered works: "Based on UG4 (www.ug4.org/license)".
 * 
 * (3) The following bibliography is recommended for citation and must be
 * preserved in all covered files:
 * "Reiter, S., Vogel, A., Heppner, I., Rupp, M., and Wittum, G. A massively
 *   parallel geometric multigrid solver on hierarchically distributed grids.
 *   Computing and visualization in science 16, 4 (2013), 151-164"
 * "Vogel, A., Reiter, S., Rupp, M., Nägel, A., and Wittum, G. UG4 -- a novel
 *   flexible software system for simulating pde based models on high performance
 *   computers. Computing and visualization in science 16, 4 (2013), 165-179"
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * GNU Lesser General Public License for more details.
 */

#ifndef __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA_IMPL_
#define __H__UG_BRIDGE__BRIDGES__USER_DATA__USER_DATA_IMPL_

#include "lua_user_data.h"


extern "C" {
#include "bindings/lua/externals/lua/lauxlib.h"
}

#include "lib_disc/spatial_disc/user_data/linker/linker_traits.h"
#include "lib_disc/spatial_disc/user_data/const_user_data.h"


#include "info_commands.h"
#include "common/util/number_util.h"
#include "bindings/lua/lua_traits.h"

//#if 0
// #define PROFILE_CALLBACK() PROFILE_FUNC_GROUP("luacallback")
//#define PROFILE_CALLBACK_BEGIN(name) PROFILE_BEGIN_GROUP(name, "luacallback")
// #define PROFILE_CALLBACK_END()	PROFILE_END()
//#else
#define PROFILE_CALLBACK()
#define PROFILE_CALLBACK_BEGIN(name)
#define PROFILE_CALLBACK_END()
//#endif

namespace ug {

#ifdef USE_LUA2C
	extern bool useLuaCompiler;
#endif



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
//	get lua state
	m_L = script::GetDefaultLuaState();

//	obtain a reference
	lua_getglobal(m_L, m_callbackName.c_str());

//	make sure that the reference is valid
	if(lua_isnil(m_L, -1)){
		UG_THROW(name() << ": Specified lua callback "
						"does not exist: " << m_callbackName);
	}

//	store reference to lua function
	m_callbackRef = luaL_ref(m_L, LUA_REGISTRYINDEX);

//	make a test run
	check_callback_returns(m_L, m_callbackRef, m_callbackName.c_str(), true);
	
	#ifdef USE_LUA2C
		if(useLuaCompiler) m_luaComp.create(luaCallback);
	#endif
}

template <typename TData, int dim, typename TRet>
LuaUserData<TData,dim,TRet>::LuaUserData(LuaFunctionHandle handle)
	: m_callbackName("__anonymous__lua__function__"), m_bFromFactory(false)
{
//	get lua state
	m_L = script::GetDefaultLuaState();

//	store reference to lua function
	m_callbackRef = handle.ref;

//	make a test run
	check_callback_returns(m_L, m_callbackRef, m_callbackName.c_str(), true);

	#ifdef USE_LUA2C
//		UG_THROW("LuaFunctionHandle usage currently not supported with LUA2C.");
		if(useLuaCompiler) m_luaComp.create(m_callbackName.c_str(), &handle);
	#endif
}


template <typename TData, int dim, typename TRet>
bool LuaUserData<TData,dim,TRet>::
check_callback_returns(lua_State* L, int callbackRef, const char* callName, const bool bThrow)
{
    PROFILE_CALLBACK()
//	get current stack level
	const int level = lua_gettop(L);

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
					"for callback\n"<<callName<< " (" << bridge::GetLUAScriptFunctionDefined(callName) << ")"
							"\nRequired: "<<retSize<<", passed: "<<numResults
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
					"for callback\n"<<callName<< " (" << bridge::GetLUAScriptFunctionDefined(callName) << ")"
							"\nUse signature as follows:\n"
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
					"for callback\n"<<callName<< " (" << bridge::GetLUAScriptFunctionDefined(callName) << ")"
							"\nUse signature as follows:\n"
							<< signature());
		}
		else{
			bRet = false;
		}
	}

//	pop values
	lua_pop(L, numResults);

//	return match
	return bRet;
}

template <typename TData, int dim, typename TRet>
bool LuaUserData<TData,dim,TRet>::
check_callback_returns(LuaFunctionHandle handle, const bool bThrow)
{
    PROFILE_CALLBACK()
//	get lua state
	lua_State* L = script::GetDefaultLuaState();

//	forward call
	bool bRet = check_callback_returns(L, handle.ref, "__lua_function_handle__", bThrow);

//	return match
	return bRet;
}

template <typename TData, int dim, typename TRet>
bool LuaUserData<TData,dim,TRet>::
check_callback_returns(const char* callName, const bool bThrow)
{
    PROFILE_CALLBACK()
//	get lua state
	lua_State* L = script::GetDefaultLuaState();

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

//	forward call
	bool bRet = check_callback_returns(L, callbackRef, callName, bThrow);

//	free reference to callback
	luaL_unref(L, LUA_REGISTRYINDEX, callbackRef);

//	return match
	return bRet;
}

template <typename TData, int dim, typename TRet>
TRet LuaUserData<TData,dim,TRet>::
evaluate(TData& D, const MathVector<dim>& x, number time, int si) const
{
    PROFILE_CALLBACK()
    #ifdef USE_LUA2C
	if(useLuaCompiler && m_luaComp.is_valid())
	{
		double d[dim+2];
		for(int i=0; i<dim; i++)
			d[i] = x[i];
		d[dim] = time;	
		d[dim+1] = si;	
		double ret[lua_traits<TData>::size+1];
		m_luaComp.call(ret, d);
		//TData D2;
		TRet *t=nullptr;
		lua_traits<TData>::read(D, ret, t);
		return lua_traits<TRet>::do_return(ret[0]);
	}
	else
	#endif
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
						"callback '" << m_callbackName << "'.\n"
						"Use signature as follows:\n"
						<< signature());

	//	pop values
		lua_pop(m_L, retSize);

	//	forward flag
		return lua_traits<TRet>::do_return(res);
	}
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
LuaUserDataFactory<TData,dim,TRet>::provide_or_create(const std::string& name)
{
	PROFILE_CALLBACK();
	using Map = std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> >;
	using iterator = typename Map::iterator;

//	check for element
	iterator iter = m_mData.find(name);

//	if name does not exist, create new one
	if(iter == m_mData.end())
	{
		SmartPtr<LuaUserData<TData,dim,TRet> > sp
			= make_sp(new LuaUserData<TData,dim,TRet>(name.c_str()));

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
LuaUserDataFactory<TData,dim,TRet>::remove(const std::string& name)
{
	using Map = std::map<std::string, std::pair<LuaUserData<TData,dim,TRet>*, int*> >;
	using iterator = typename Map::iterator;

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
	m_L = script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
	m_cbDerivRef.clear();
	m_cbDerivName.clear();
	set_lua_value_callback(luaCallback, numArgs);
	#ifdef USE_LUA2C
		if(useLuaCompiler) m_luaComp.create(luaCallback);
	#endif
}

template <typename TData, int dim, typename TDataIn>
LuaUserFunction<TData,dim,TDataIn>::
LuaUserFunction(const char* luaCallback, size_t numArgs, bool bPosTimeNeed)
	: m_numArgs(numArgs), m_bPosTimeNeed(bPosTimeNeed)
{
	m_L = script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
	m_cbDerivRef.clear();
	m_cbDerivName.clear();
	set_lua_value_callback(luaCallback, numArgs);
	#ifdef USE_LUA2C
		m_luaComp_Deriv.clear();
	#endif
}


template <typename TData, int dim, typename TDataIn>
LuaUserFunction<TData,dim,TDataIn>::
LuaUserFunction(LuaFunctionHandle handle, size_t numArgs)
	: m_numArgs(numArgs), m_bPosTimeNeed(false)
{
	m_L = script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
	m_cbDerivRef.clear();
	m_cbDerivName.clear();
	set_lua_value_callback(handle, numArgs);
	#ifdef USE_LUA2C
		if(useLuaCompiler){
			UG_LOG("WARNING (in LuaUserFunction): LUA2C compiler "
					"can't be executed for FunctionHandle.\n");
		}
	#endif
}

template <typename TData, int dim, typename TDataIn>
LuaUserFunction<TData,dim,TDataIn>::
LuaUserFunction(LuaFunctionHandle handle, size_t numArgs, bool bPosTimeNeed)
	: m_numArgs(numArgs), m_bPosTimeNeed(bPosTimeNeed)
{
	m_L = script::GetDefaultLuaState();
	m_cbValueRef = LUA_NOREF;
	m_cbDerivRef.clear();
	m_cbDerivName.clear();
	set_lua_value_callback(handle, numArgs);
	#ifdef USE_LUA2C
		m_luaComp_Deriv.clear();
	#endif
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

	#ifdef USE_LUA2C
		m_luaComp_Deriv.resize(numArgs);
	#endif
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::
set_lua_value_callback(LuaFunctionHandle handle, size_t numArgs)
{
//	store name (string) of callback
	m_cbValueName = "__anonymous__lua__function__";

//	if a callback was already set, we have to free the old one
	free_callback_ref();

//	store reference to lua function
	m_cbValueRef = handle.ref;

//	remember number of arguments to be used
	m_numArgs = numArgs;
	m_cbDerivName.resize(numArgs);
	m_cbDerivRef.resize(numArgs, LUA_NOREF);

//	set num inputs for linker
	set_num_input(numArgs);

	#ifdef USE_LUA2C
		m_luaComp_Deriv.resize(numArgs);
	#endif
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
	
	#ifdef USE_LUA2C
		if(useLuaCompiler) m_luaComp_Deriv[arg].create(luaCallback);
	#endif
	
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::set_deriv(size_t arg, LuaFunctionHandle handle)
{
//	check number of arg
	if(arg >= m_numArgs)
		UG_THROW("LuaUserFunction::set_lua_deriv_callback: Trying "
				"to set a derivative for argument " << arg <<", that "
				"does not exist. Number of arguments is "<<m_numArgs);

//	store name (string) of callback
	m_cbDerivName[arg] = std::string("__anonymous__lua__function__");

//	free old reference
	free_deriv_callback_ref(arg);

//	store reference to lua function
	m_cbDerivRef[arg] = handle.ref;

	#ifdef USE_LUA2C
	//	if(useLuaCompiler) m_luaComp_Deriv[arg].create(luaCallback);
	#endif

}




template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::operator () (TData& out, int numArgs, ...) const
{
    PROFILE_CALLBACK();
    #ifdef USE_LUA2C
	if(useLuaCompiler && m_luaComp.is_valid())
	{
		double d[20];
		//	get list of arguments
		va_list ap2;
		va_start(ap2, numArgs);

	//	read all arguments and push them to the lua stack
		for(int i = 0; i < numArgs; ++i)
			d[i] = va_arg(ap2, double);
		va_end(ap2);

		double ret[lua_traits<TData>::size+1];

		UG_ASSERT(m_luaComp.num_in() == numArgs && m_luaComp.num_out() == lua_traits<TData>::size,
				m_luaComp.name() << ", " << m_luaComp.num_in() << " != " << numArgs << " or " << m_luaComp.num_out() << " != " << lua_traits<TData>::size);
		m_luaComp.call(ret, d);
		//TData D2;
		void *t=nullptr;
		//TData out2;
		lua_traits<TData>::read(out, ret, t);
		return;
	}
	else
	#endif
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
			UG_COND_THROW(IsFiniteAndNotTooBig(out)==false, out);
		}
		UG_CATCH_THROW("LuaUserFunction::operator(...): Error while running "
						"callback '" << m_cbValueName << "'");

	//	pop values
		lua_pop(m_L, retSize);
	}
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::eval_value(TData& out, const std::vector<TDataIn>& dataIn,
													const MathVector<dim>& x, number time, int si) const
{
    PROFILE_CALLBACK();
    #ifdef USE_LUA2C
	if(useLuaCompiler && m_luaComp.is_valid())
	{
		double d[20];

	//	read all arguments and push them to the lua stack
		for(size_t i = 0; i < dataIn.size(); ++i)
			d[i] = dataIn[i];
		if(m_bPosTimeNeed){
			for(int i=0; i<dim; i++)
				d[i+m_numArgs] = x[i];
			d[dim+m_numArgs]=time;
			d[dim+m_numArgs+1]=si;
			UG_ASSERT(dim+m_numArgs+1 < 20, m_luaComp.name());
		}

		double ret[lua_traits<TData>::size];
		m_luaComp.call(ret, d);
		//TData D2;
		void *t=nullptr;
		//TData out2;
		UG_ASSERT(m_luaComp.num_out() == lua_traits<TData>::size, m_luaComp.name() << ", " << m_luaComp.num_out() << " != " << lua_traits<TData>::size);
		lua_traits<TData>::read(out, ret, t);
		return;
	}
	else
	#endif
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
			UG_COND_THROW(IsFiniteAndNotTooBig(out)==false, out);
		}
		UG_CATCH_THROW("LuaUserFunction::eval_value(...): Error while "
						"running callback '" << m_cbValueName << "'");

	//	pop values
		lua_pop(m_L, retSize);
	}
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::eval_deriv(TData& out, const std::vector<TDataIn>& dataIn,
		 	 	 	 	 	 	 	 	 	 	 	 const MathVector<dim>& x, number time, int si, size_t arg) const
{
	PROFILE_CALLBACK();
	#ifdef USE_LUA2C
	if(useLuaCompiler && m_luaComp_Deriv[arg].is_valid()
		&& dim+m_numArgs+1 < 20 && m_luaComp_Deriv[arg].num_out() == lua_traits<TData>::size)
	{
		const bridge::LUACompiler &luaComp = m_luaComp_Deriv[arg];
		double d[25];
		UG_ASSERT(dim+m_numArgs+1 < 20, luaComp.name());
		for(size_t i=0; i<m_numArgs; i++)
			d[i] = dataIn[i];
		if(m_bPosTimeNeed){
			for(int i=0; i<dim; i++)
				d[i+m_numArgs] = x[i];
			d[dim+m_numArgs]=time;
			d[dim+m_numArgs+1]=si;
			UG_ASSERT(dim+m_numArgs+1 < 20, luaComp.name());
		}
		UG_ASSERT(luaComp.num_out() == lua_traits<TData>::size,
				luaComp.name() << " has wrong number of outputs: is " << luaComp.num_out() << ", needs " << lua_traits<TData>::size);
		double ret[lua_traits<TData>::size+1];
		luaComp.call(ret, d);
		//TData D2;
		void *t=nullptr;
		//TData out2;
		lua_traits<TData>::read(out, ret, t);
		return;
	}
	else
	#endif
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
			UG_COND_THROW(IsFiniteAndNotTooBig(out)==false, out);
		}
		UG_CATCH_THROW("LuaUserFunction::eval_deriv(...): Error while "
					"running callback '" << m_cbDerivName[arg] << "'");

	//	pop values
		lua_pop(m_L, retSize);
	}
}


template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::
evaluate (TData& value,
          const MathVector<dim>& globIP,
          number time, int si) const
{
	PROFILE_CALLBACK();
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t c = 0; c < vDataIn.size(); ++c)
		(*m_vpUserData[c])(vDataIn[c], globIP, time, si);

//	evaluate data at ip
	eval_value(value, vDataIn, globIP, time, si);

	UG_COND_THROW(IsFiniteAndNotTooBig(value)==false, value);
}

template <typename TData, int dim, typename TDataIn>
template <int refDim>
void LuaUserFunction<TData,dim,TDataIn>::
evaluate(TData vValue[],
         const MathVector<dim> vGlobIP[],
         number time, int si,
         GridObject* elem,
         const MathVector<dim> vCornerCoords[],
         const MathVector<refDim> vLocIP[],
         const size_t nip,
         LocalVector* u,
         const MathMatrix<refDim, dim>* vJT) const
{
	PROFILE_CALLBACK();
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

//	gather all input data for this ip
	for(size_t ip = 0; ip < nip; ++ip)
	{
		for(size_t c = 0; c < vDataIn.size(); ++c)
			(*m_vpUserData[c])(vDataIn[c], vGlobIP[ip], time, si, elem, vCornerCoords, vLocIP[ip], u);

	//	evaluate data at ip
		eval_value(vValue[ip], vDataIn, vGlobIP[ip], time, si);
		UG_COND_THROW(IsFiniteAndNotTooBig(vValue[ip])==false, vValue[ip]);
	}
}

template <typename TData, int dim, typename TDataIn>
template <int refDim>
void LuaUserFunction<TData,dim,TDataIn>::
eval_and_deriv(TData vValue[],
				const MathVector<dim> vGlobIP[],
				number time, int si,
				GridObject* elem,
				const MathVector<dim> vCornerCoords[],
				const MathVector<refDim> vLocIP[],
				const size_t nip,
				LocalVector* u,
				bool bDeriv,
				int s,
				std::vector<std::vector<TData> > vvvDeriv[],
				const MathMatrix<refDim, dim>* vJT)
{
	PROFILE_CALLBACK();
//	vector of data for all inputs
	std::vector<TDataIn> vDataIn(this->num_input());

	for(size_t ip = 0; ip < nip; ++ip)
	{
	//	gather all input data for this ip
		for(size_t c = 0; c < vDataIn.size(); ++c)
			vDataIn[c] = m_vpUserData[c]->value(this->series_id(c,s), ip);

	//	evaluate data at ip
		eval_value(vValue[ip], vDataIn, vGlobIP[ip], time, si);
	}

//	check if derivative is required
	if(!bDeriv || this->zero_derivative()) return;

//	clear all derivative values
	this->set_zero(vvvDeriv, nip);

//	loop all inputs
	for(size_t c = 0; c < vDataIn.size(); ++c)
	{
	//	check if we have the derivative w.r.t. this input, and the input has derivative
		if(m_cbDerivRef[c] == LUA_NOREF || m_vpUserData[c]->zero_derivative()) continue;

	//	loop ips
		for(size_t ip = 0; ip < nip; ++ip)
		{
		//	gather all input data for this ip
			for(size_t i = 0; i < vDataIn.size(); ++i)
				vDataIn[i] = m_vpUserData[i]->value(this->series_id(c,s), ip); //< series_id(c,s) or series_id(i,s)

		//	data of derivative w.r.t. one component at ip-values
			TData derivVal;

		//	evaluate data at ip
			eval_deriv(derivVal, vDataIn, vGlobIP[ip], time, si, c);

		//	loop functions
			for(size_t fct = 0; fct < this->input_num_fct(c); ++fct)
			{
			//	get common fct id for this function
				const size_t commonFct = this->input_common_fct(c, fct);

			//	loop dofs
				for(size_t dof = 0; dof < this->num_sh(fct); ++dof)
				{
					linker_traits<TData, TDataIn>::
					mult_add(vvvDeriv[ip][commonFct][dof],
							 derivVal,
							 m_vpDependData[c]->deriv(this->series_id(c,s), ip, fct, dof));
					UG_COND_THROW(IsFiniteAndNotTooBig(vvvDeriv[ip][commonFct][dof])==false, vvvDeriv[ip][commonFct][dof]);
				}
			}
		}
	}
}

/**
 * TODO: Note this is a public (non-virtual) function whose argument
 * should be consistent with the number of the arguments. Should not it also
 * resize the array for the references to the derivatives?
 */
template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::set_num_input(size_t num)
{
//	resize arrays
	m_vpUserData.resize(num);
	m_vpDependData.resize(num);

//	forward size to base class
	base_type::set_num_input(num);
}

template <typename TData, int dim, typename TDataIn>
void LuaUserFunction<TData,dim,TDataIn>::
set_input(size_t i, SmartPtr<CplUserData<TDataIn, dim> > data)
{
	UG_ASSERT(i < m_vpUserData.size(), "Input not needed");
	UG_ASSERT(i < m_vpDependData.size(), "Input not needed");

//	check input number
	if(i >= this->num_input())
		UG_THROW("LuaUserFunction::set_input: Only " << this->num_input()
						<< " inputs can be set. Use 'set_num_input' to increase"
						" the number of needed inputs.");

//	remember userdata
	m_vpUserData[i] = data;

//	cast to dependent data
	m_vpDependData[i] = data.template cast_dynamic<DependentUserData<TDataIn, dim> >();

//	forward to base class
	base_type::set_input(i, data, data);
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
	m_L = script::GetDefaultLuaState();
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
void LuaFunction<TData,TDataIn>::operator () (TData& out, int numArgs, ...)
{
	PROFILE_CALLBACK_BEGIN(operatorBracket);
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
		UG_COND_THROW(IsFiniteAndNotTooBig(out)==false, out);
	}
	UG_CATCH_THROW("LuaFunction::operator(...): Error while running "
					"callback '" << m_cbValueName << "'");

//	pop values
	lua_pop(m_L, retSize);

    PROFILE_CALLBACK_END();
}



} // end namespace ug

#endif
